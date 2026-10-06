//! Pure Rust implementation of [Sparc](https://github.com/keithyin/Sparc), a
//! sparsity-based consensus algorithm for long erroneous sequencing reads
//! (Ye C, Ma Z. 2016, PeerJ 4:e2016).
//!
//! All interaction goes through [`sparc_consensus`]: give it a backbone
//! sequence plus alignments ([`Query`]) of reads against that backbone, and it
//! returns the polished consensus. Alignments can be built by hand or parsed
//! from blasr m5 rows/files ([`Query::from_m5_row`], [`parse_m5`],
//! [`M5Reader`]).
//!
//! # Example
//!
//! ```
//! use sparc::{Query, SparcConfig, sparc_consensus};
//!
//! let backbone = "GATCGGGCTAA";
//! let queries = ["GATCGCGCTAA", "GATCGCGCCAA", "GATCGCGCCAA", "GATCGCGCCAA"]
//!     .iter()
//!     .map(|seq| Query::new(backbone, *seq, 0, backbone.len()))
//!     .collect::<Vec<_>>();
//!
//! let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
//! assert!(!consensus.seq.is_empty());
//! ```
//!
//! Inputs are validated before the algorithm runs; invalid input returns
//! [`SparcError`] instead of misbehaving.

#![warn(missing_docs)]

mod best_path;
mod config;
mod encoding;
mod error;
mod graph;
mod m5;
mod normalize;
mod pipeline;
mod validate;

pub use config::{DebugConfig, ScoringMethod, SparcConfig};
pub use error::SparcError;
pub use m5::{M5Reader, parse_m5};

use validate::validate_inputs;

/// 最优路径在 backbone 上的区间语义（[`Consensus::region`] 的返回值）。
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PathRegion {
    /// 路径首尾都是 backbone 节点：节点下标区间 `[start, end)`。
    Backbone {
        /// 区间起点下标（含）。
        start: usize,
        /// 区间终点下标（不含）。
        end: usize,
    },
    /// 路径头部不在 backbone 上（例如起点是插入分支/孤儿链首节点），
    /// 仅知终点下标 `end`（含）。
    OrphanHead {
        /// 路径终点下标（含）。
        end: usize,
    },
    /// 无可信路径（例如无 query 输入）：`seq` 为原 backbone。
    Fallback,
}

/// consensus 结果。
///
/// `start`/`end` 是最优路径在 backbone 上的节点下标区间 `[start, end)`：
/// - `start == None`：路径头部不在 backbone 上（例如起点是一个插入分支节点）；
/// - `end == 0`：无可信路径（如无 query 输入），此时 `seq` 为原 backbone。
///
/// 这组哨兵编码可用 [`Consensus::region`] 以 [`PathRegion`] 形式直接获取。
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Consensus {
    /// consensus 序列（fallback 时为原 backbone）。
    pub seq: String,
    /// 最优路径起点下标；`None` 表示路径头部不在 backbone 上（孤儿链首）。
    pub start: Option<usize>,
    /// 最优路径终点下标（不含）；`0` 表示无可信路径，`seq` 原样返回 backbone。
    pub end: usize,
}

impl Consensus {
    /// 最优路径区间的语义化视图（与 `start`/`end` 哨兵编码一一对应）。
    pub fn region(&self) -> PathRegion {
        match (self.start, self.end) {
            (Some(start), end) => PathRegion::Backbone { start, end },
            (None, 0) => PathRegion::Fallback,
            (None, end) => PathRegion::OrphanHead { end },
        }
    }
}

/// 一条 read 相对 backbone 的比对，对应一行 blasr m5 记录。
///
/// `query_aligned_sequence` / `target_aligned_sequence` 是带 `-` gap 列的比对串，
/// 两者必须等长；`target_aligned_sequence` 中非 `-` 的碱基数必须等于
/// `target_end - target_start`。以上约束由 [`sparc_consensus`] 在调用时校验。
///
/// # Example
///
/// ```
/// use sparc::Query;
///
/// // 参数顺序：先 target（backbone 方向）后 query；坐标为 backbone 上的 [start, end)
/// let query = Query::new("GATC-GGCTAA", "GATCGCGCTAA", 0, 11).reverse_strand();
/// assert!(query.is_reverse_strand());
/// assert_eq!(query.target_start(), 0);
/// assert_eq!(query.target_aligned_sequence(), "GATC-GGCTAA");
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Query {
    pub(crate) query_aligned_sequence: String,
    pub(crate) target_aligned_sequence: String,
    pub(crate) is_reverse_strand: bool,
    /// 来自 m5 的 qStart。仅信息性，算法不使用。
    pub(crate) query_start: usize,
    /// 来自 m5 的 qEnd。仅信息性，算法不使用。
    pub(crate) query_end: usize,
    pub(crate) target_start: usize,
    pub(crate) target_end: usize,
}

impl Query {
    /// 由 target/query 两条比对串构造（参数顺序：先 target 后 query）。
    /// `target_aligned` 为 backbone 方向的比对串，坐标为 backbone 上的
    /// `[target_start, target_end)`。默认按正链（tStrand = '+'）处理。
    pub fn new(
        target_aligned: impl Into<String>,
        query_aligned: impl Into<String>,
        target_start: usize,
        target_end: usize,
    ) -> Self {
        Self {
            query_aligned_sequence: query_aligned.into(),
            target_aligned_sequence: target_aligned.into(),
            is_reverse_strand: false,
            query_start: 0,
            query_end: 0,
            target_start,
            target_end,
        }
    }

    /// 标记该比对来自负链（m5 的 tStrand == '-'）。
    ///
    /// 传入的两条比对串保持 m5 负链行的原始方向（read 方向），
    /// 算法内部会对其做反向互补。
    #[must_use = "reverse_strand 返回标记后的新 Query，丢弃返回值则标记无效"]
    pub fn reverse_strand(mut self) -> Self {
        self.is_reverse_strand = true;
        self
    }

    /// 解析一行 blasr m5 记录（19 个空白分隔字段，见 README）。
    ///
    /// 只取绑定层用到的字段：qStart/qEnd、tStart/tEnd、tStrand、
    /// qAlignedSeq、tAlignedSeq（matchPattern 等忽略）。
    /// 字母表/坐标一致性由 [`sparc_consensus`] 校验。
    pub fn from_m5_row(row: &str) -> Result<Self, SparcError> {
        m5::parse_m5_row(row)
    }

    /// read 方向的比对串（带 `-` gap 列）。
    pub fn query_aligned_sequence(&self) -> &str {
        &self.query_aligned_sequence
    }

    /// backbone 方向的比对串（带 `-` gap 列）。
    pub fn target_aligned_sequence(&self) -> &str {
        &self.target_aligned_sequence
    }

    /// 是否为负链（m5 tStrand == '-'）。
    pub fn is_reverse_strand(&self) -> bool {
        self.is_reverse_strand
    }

    /// read 上比对起点的 m5 坐标（qStart）。仅信息性，算法不使用。
    pub fn query_start(&self) -> usize {
        self.query_start
    }

    /// read 上比对终点的 m5 坐标（qEnd）。仅信息性，算法不使用。
    pub fn query_end(&self) -> usize {
        self.query_end
    }

    /// backbone 上比对区间起点（含）。
    pub fn target_start(&self) -> usize {
        self.target_start
    }

    /// backbone 上比对区间终点（不含）。
    pub fn target_end(&self) -> usize {
        self.target_end
    }
}

/// 由 backbone 和一组比对计算 consensus。
///
/// 返回 [`Consensus`]，其 `start`/`end` 为最优路径在 backbone 上的
/// `[start, end)` 节点下标区间（`end == 0` 表示回退为原 backbone，
/// `start == None` 表示路径头部不在 backbone 上），也可用
/// [`Consensus::region`] 获取语义化视图。
///
/// 所有输入在进入算法前完成校验；不满足约束时返回 [`SparcError`]。
pub fn sparc_consensus(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<Consensus, SparcError> {
    validate_inputs(backbone, queries, config)?;
    pipeline::compute_consensus(backbone, queries, config)
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    fn make_queries(backbone: &str) -> Vec<Query> {
        let q = |query_aligned_sequence: &str| {
            Query::new(backbone, query_aligned_sequence, 0, backbone.len())
        };
        vec![
            q("GATCGCGCTAA"),
            q("GATCGCGCCAA"),
            q("GCTCGGCCCAA"),
            q("GCTCGGCCCAA"),
            q("GCTCGGCCCAA"),
            q("GATCGCGCCAA"),
            q("GATCGCGCCAA"),
        ]
    }

    fn reverse_complement(sequence: &str) -> String {
        sequence
            .bytes()
            .rev()
            .map(|base| match base {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                b'G' => 'C',
                _ => '-',
            })
            .collect()
    }

    #[test]
    fn test_sparc_consensus() {
        let backbone = "GATCGGGCTAA";
        let queries = make_queries(backbone);

        let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
        assert!(!consensus.seq.is_empty());
        assert_eq!(consensus.start, Some(0));
        assert!(consensus.end > 0);
    }

    /// 回归：空 queries 时必须走 fallback 路径，原样返回 backbone
    /// （seq == backbone、start == None、end == 0、region == Fallback）。
    #[test]
    fn test_empty_queries_returns_backbone() {
        let backbone = "GATCGGGCTAA";

        let consensus = sparc_consensus(backbone, &[], &SparcConfig::default()).unwrap();
        assert_eq!(consensus.seq, backbone);
        assert_eq!(consensus.start, None);
        assert_eq!(consensus.end, 0);
        assert_eq!(consensus.region(), PathRegion::Fallback);
    }

    /// 回归：重复调用 50 次结果稳定，不存在跨调用残留状态。
    #[test]
    fn test_sparc_consensus_repeated() {
        let backbone = "GATCGGGCTAA";
        let queries = make_queries(backbone);

        for _ in 0..50 {
            let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
            assert!(!consensus.seq.is_empty());
            assert!(consensus.start.expect("start should be set") < consensus.end);
        }
    }

    /// 回归：3' 端插入使非 backbone 节点挂在最后一个 backbone 节点上，
    /// 该分支路径必须参与最优路径计算且不破坏结果。
    #[test]
    fn test_consensus_with_terminal_insertion() {
        let backbone = "GATCGGGCTAA";
        let mut queries = make_queries(backbone);

        // target 比对串末尾为 '-'，即在 backbone 末端插入一个碱基
        queries.push(Query::new(
            format!("{backbone}-"),
            "GATCGGGCTAAA",
            0,
            backbone.len(),
        ));
        queries.push(Query::new(
            format!("{backbone}-"),
            "GATCGGGCTAAG",
            0,
            backbone.len(),
        ));

        let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
        assert!(!consensus.seq.is_empty());
        assert!(consensus.start.expect("start should be set") < consensus.end);
    }

    /// 回归：比对起始处即 mismatch 的 read 会在链首产生孤儿节点
    /// （不挂在任何 backbone 节点的右子图内），结果仍须非空。
    #[test]
    fn test_leading_mismatch_queries() {
        let backbone = "GATCGGGCTAA";
        let queries: Vec<Query> = (0..1000)
            .map(|_| Query::new(backbone, "TATCGGGCTAA", 0, backbone.len()))
            .collect();

        let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
        assert!(!consensus.seq.is_empty());
    }

    /// 负链：m5 负链行的两条比对串按 read 方向给出，算法内部反向互补后
    /// 应与直接给正向比对得到完全一致的结果。
    #[test]
    fn test_reverse_strand_matches_forward() {
        let backbone = "GATCGGGCTAA";

        let forward = make_queries(backbone);
        let reversed: Vec<Query> = forward
            .iter()
            .map(|q| {
                Query::new(
                    reverse_complement(q.target_aligned_sequence()),
                    reverse_complement(q.query_aligned_sequence()),
                    q.target_start(),
                    q.target_end(),
                )
                .reverse_strand()
            })
            .collect();

        let a = sparc_consensus(backbone, &forward, &SparcConfig::default()).unwrap();
        let b = sparc_consensus(backbone, &reversed, &SparcConfig::default()).unwrap();
        assert_eq!(a, b);
    }

    /// threshold 行为回归：自适应阈值必须真正影响结果。
    /// threshold=100 使所有边的 new_score 为负、从不更新，全部节点 score=0，
    /// 触发 fallback（原样返回 backbone）；threshold=-0.1 则正常产出 consensus。
    #[test]
    fn test_threshold_changes_behavior() {
        let backbone = "GATCGGGCTAA";
        let queries = make_queries(backbone);

        let config = SparcConfig {
            threshold: -0.1,
            ..SparcConfig::default()
        };
        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert_ne!(consensus.seq, backbone);
        assert!(consensus.end > 0);

        let config = SparcConfig {
            threshold: 100.0,
            ..SparcConfig::default()
        };
        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert_eq!(consensus.seq, backbone);
        assert_eq!(consensus.end, 0);
        assert_eq!(consensus.start, None);
    }

    /// 回归：k=3 且 cov_radius 大于 backbone 节点数时，滑动窗口半径
    /// 必须被正确 clamp（历史上此处越界读 node_vec），输出确定性。
    #[test]
    fn test_k3_large_radius() {
        let bb = "GATCGGGCTAAACGTACGATCGATCGATCGGATCCGATACGTACGTACGATCGTACGATC";
        let config = SparcConfig {
            kmer: 3,
            cov_radius: 500,
            ..SparcConfig::default()
        };
        let exact = Query::new(bb, bb, 0, bb.len());
        let mut queries = vec![exact.clone(), exact.clone(), exact];
        // 首列替换：产生孤儿链首；末端插入：产生末端分支
        queries.push(Query::new(bb, format!("T{}", &bb[1..]), 0, bb.len()));
        queries.push(Query::new(format!("{bb}-"), format!("{bb}A"), 0, bb.len()));

        let first = sparc_consensus(bb, &queries, &config).unwrap();
        assert!(!first.seq.is_empty());
        // 确定性：同输入重复调用结果一致
        for _ in 0..3 {
            let again = sparc_consensus(bb, &queries, &config).unwrap();
            assert_eq!(again, first);
        }
    }

    /// `Consensus::region` 与 `start`/`end` 哨兵编码的三态映射。
    #[test]
    fn test_consensus_region() {
        assert_eq!(
            Consensus {
                seq: String::new(),
                start: Some(2),
                end: 5
            }
            .region(),
            PathRegion::Backbone { start: 2, end: 5 }
        );
        assert_eq!(
            Consensus {
                seq: String::new(),
                start: None,
                end: 3
            }
            .region(),
            PathRegion::OrphanHead { end: 3 }
        );
        assert_eq!(
            Consensus {
                seq: String::new(),
                start: None,
                end: 0
            }
            .region(),
            PathRegion::Fallback
        );

        // 端到端：正常路径 → Backbone；空 queries → Fallback
        let backbone = "GATCGGGCTAA";
        let normal =
            sparc_consensus(backbone, &make_queries(backbone), &SparcConfig::default()).unwrap();
        assert_eq!(
            normal.region(),
            PathRegion::Backbone {
                start: normal.start.expect("start should be set"),
                end: normal.end
            }
        );
    }

    #[test]
    fn test_validation_errors() {
        let backbone = "GATCGGGCTAA";
        let config = SparcConfig::default();
        let q = |query_aligned: &str, target_aligned: &str, span_start: usize, span_end: usize| {
            Query::new(target_aligned, query_aligned, span_start, span_end)
        };

        // kmer 越界（下界 / 上界）
        for kmer in [0, 17] {
            let cfg = SparcConfig {
                kmer,
                ..SparcConfig::default()
            };
            assert_eq!(
                sparc_consensus(backbone, &[], &cfg),
                Err(SparcError::InvalidKmer { kmer })
            );
        }

        // backbone 长度不足
        let cfg = SparcConfig {
            kmer: 5,
            ..SparcConfig::default()
        };
        assert_eq!(
            sparc_consensus("GAT", &[], &cfg),
            Err(SparcError::BackboneTooShort {
                backbone_len: 3,
                kmer: 5
            })
        );

        // backbone 含 N
        assert_eq!(
            sparc_consensus("GATCGNGCTAA", &[], &config),
            Err(SparcError::InvalidBackboneBase {
                base: 'N',
                position: 5
            })
        );

        // threshold 为 NaN / +inf（NaN != NaN，用 matches! 而非 assert_eq! 比较）
        let cfg = SparcConfig {
            threshold: f64::NAN,
            ..SparcConfig::default()
        };
        assert!(matches!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidThreshold { .. })
        ));
        let cfg = SparcConfig {
            threshold: f64::INFINITY,
            ..SparcConfig::default()
        };
        assert!(matches!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidThreshold { .. })
        ));

        // threshold 为 -inf 合法：与其他负值同义（关闭自适应）
        let cfg = SparcConfig {
            threshold: f64::NEG_INFINITY,
            ..SparcConfig::default()
        };
        assert!(sparc_consensus(backbone, &make_queries(backbone), &cfg).is_ok());

        // 空比对
        assert_eq!(
            sparc_consensus(backbone, &[q("", "", 0, 11)], &config),
            Err(SparcError::EmptyAlignment)
        );

        // 长度不齐
        assert_eq!(
            sparc_consensus(backbone, &[q("GATC", "GATCG", 0, 5)], &config),
            Err(SparcError::AlignedLengthMismatch {
                query_len: 4,
                target_len: 5
            })
        );

        // 非法碱基
        assert_eq!(
            sparc_consensus(backbone, &[q("GATCN", "GATCN", 0, 5)], &config),
            Err(SparcError::InvalidAlignedBase {
                base: 'N',
                position: 4
            })
        );

        // span 越界
        assert_eq!(
            sparc_consensus(backbone, &[q("GATCG", "GATCG", 0, 99)], &config),
            Err(SparcError::InvalidTargetSpan {
                target_start: 0,
                target_end: 99,
                backbone_len: backbone.len(),
            })
        );

        // span 与非 gap 碱基数不一致（含一个 gap 列）
        assert_eq!(
            sparc_consensus(backbone, &[q("GATC-GGCTAA", "GATC-GGCTAA", 0, 11)], &config),
            Err(SparcError::SpanMismatch {
                span: 11,
                target_bases: 10
            })
        );
    }

    #[test]
    fn test_parse_m5_row() {
        let row = "read_24/0_10 10 0 12 + ref_template 11 0 11 + 100 10 0 0 0 60 \
                   GATCGCGCTAA 10M GATCGGGCTAA";
        let q = Query::from_m5_row(row).unwrap();
        assert_eq!(q.query_aligned_sequence(), "GATCGCGCTAA");
        assert_eq!(q.target_aligned_sequence(), "GATCGGGCTAA");
        assert_eq!(q.query_start(), 0);
        assert_eq!(q.query_end(), 12);
        assert_eq!(q.target_start(), 0);
        assert_eq!(q.target_end(), 11);
        assert!(!q.is_reverse_strand());

        let negative_row = "read_24/0_10 10 0 12 + ref_template 11 0 11 - 100 10 0 0 0 60 \
                            GATCGCGCTAA 10M GATCGGGCTAA";
        let q = Query::from_m5_row(negative_row).unwrap();
        assert!(q.is_reverse_strand());

        let bad = Query::from_m5_row("only three fields");
        assert_eq!(
            bad.unwrap_err(),
            SparcError::M5Parse {
                line: 1,
                message: "expected 19 whitespace-separated fields, got 3".to_string(),
            }
        );

        let bad_strand = Query::from_m5_row(
            "r 10 0 12 + t 11 0 11 ? 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA",
        );
        assert!(matches!(bad_strand, Err(SparcError::M5Parse { .. })));
    }

    /// M5Reader 跳过空行且错误带 1 起始的物理行号。
    #[test]
    fn test_m5_reader_line_numbers_and_blank_lines() {
        let valid = "r 10 0 12 + t 11 0 11 + 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA";

        let data = format!("\n{valid}\n\n{valid}\n");
        let queries: Vec<Query> = M5Reader::new(Cursor::new(data))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(queries.len(), 2);
        assert_eq!(queries[0], queries[1]);

        // 行 1、3 有效，行 4 非法 → 错误定位到第 4 行
        let data = format!("{valid}\n\n{valid}\nbad line\n");
        let result: Result<Vec<Query>, SparcError> = M5Reader::new(Cursor::new(data)).collect();
        assert_eq!(
            result.unwrap_err(),
            SparcError::M5Parse {
                line: 4,
                message: "expected 19 whitespace-separated fields, got 2".to_string(),
            }
        );
    }

    #[test]
    fn test_parse_m5_testdata() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/testdata2/backbone-0.mapped.m5"
        );
        let queries =
            parse_m5(std::io::BufReader::new(std::fs::File::open(path).unwrap())).unwrap();
        assert_eq!(queries.len(), 7);
        assert!(
            queries
                .iter()
                .all(|q| q.target_start() == 0 && !q.is_reverse_strand())
        );
    }

    /// 端到端：testdata2 的 backbone + m5 应复现原版 CLI 的 out.consensus.fasta。
    #[test]
    fn test_e2e_testdata2() {
        let dir = env!("CARGO_MANIFEST_DIR");
        let fasta = std::fs::read_to_string(format!("{dir}/testdata2/backbone-0.fasta")).unwrap();
        let backbone: String = fasta.lines().skip(1).collect();
        assert_eq!(backbone, "GATCGGGCTAA");

        let f = std::fs::File::open(format!("{dir}/testdata2/backbone-0.mapped.m5")).unwrap();
        let queries = parse_m5(std::io::BufReader::new(f)).unwrap();

        // 对齐原 CLI 生成 out.consensus.fasta 时的参数：
        // k=1（默认），threshold=-0.1 关闭自适应，滑动窗口半径 200
        let config = SparcConfig {
            threshold: -0.1,
            cov_radius: 200,
            ..SparcConfig::default()
        };

        let consensus = sparc_consensus(&backbone, &queries, &config).unwrap();
        let expected =
            std::fs::read_to_string(format!("{dir}/testdata2/out.consensus.fasta")).unwrap();
        let expected_seq: String = expected.lines().skip(1).collect();
        assert_eq!(consensus.seq, expected_seq);
    }

    /// debug 输出写入指定目录（不再污染进程 CWD）。
    #[test]
    fn test_debug_output_dir() {
        let dir = tempfile::tempdir().unwrap();
        let backbone = "GATCGGGCTAA";
        let config = SparcConfig {
            debug_output: Some(DebugConfig {
                output_dir: dir.path().to_path_buf(),
                subgraph_begin: 0,
                subgraph_end: 100,
            }),
            ..SparcConfig::default()
        };

        let consensus = sparc_consensus(backbone, &make_queries(backbone), &config).unwrap();
        assert_ne!(consensus.end, 0); // 非 fallback 才有 DEBUG.consensus.fasta / subgraph_cns.dot
        for filename in [
            "align_profile.txt",
            "DEBUG.consensus.fasta",
            "subgraph.dot",
            "subgraph_cns.dot",
        ] {
            assert!(
                dir.path().join(filename).exists(),
                "{filename} not written to output_dir"
            );
        }
    }

    /// debug 输出目录不存在时，文件创建失败必须以 DebugWrite 传播而非 panic。
    #[test]
    fn test_debug_write_error_propagates() {
        let dir = tempfile::tempdir().unwrap();
        // 把一个已存在的普通文件当作输出目录，File::create 必然失败
        let blocker = dir.path().join("not-a-dir");
        std::fs::write(&blocker, b"x").unwrap();

        let config = SparcConfig {
            debug_output: Some(DebugConfig {
                output_dir: blocker,
                ..DebugConfig::default()
            }),
            ..SparcConfig::default()
        };
        assert!(matches!(
            sparc_consensus("GATCGGGCTAA", &[], &config),
            Err(SparcError::DebugWrite(_))
        ));
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_serde_roundtrip() {
        let config = SparcConfig {
            debug_output: Some(DebugConfig {
                output_dir: "/tmp/sparc-debug".into(),
                subgraph_begin: 3,
                subgraph_end: 40,
            }),
            kmer: 3,
            coverage_threshold: 5,
            scoring_method: ScoringMethod::LogRatio,
            cov_radius: 200,
            threshold: -0.1,
        };
        let json = serde_json::to_string(&config).unwrap();
        assert_eq!(serde_json::from_str::<SparcConfig>(&json).unwrap(), config);
    }
}
