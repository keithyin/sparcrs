//! Pure Rust implementation of [Sparc](https://github.com/keithyin/Sparc), a
//! sparsity-based consensus algorithm for long erroneous sequencing reads
//! (Ye C, Ma Z. 2016, PeerJ 4:e2016).
//!
//! All interaction goes through [`sparc_consensus`]: give it a backbone
//! sequence plus alignments ([`Query`]) of reads against that backbone, and it
//! returns the polished consensus.
//!
//! # Example
//!
//! ```
//! use sparc::{Query, SparcConfig, sparc_consensus};
//!
//! let backbone = "GATCGGGCTAA";
//! let queries = ["GATCGCGCTAA", "GATCGCGCCAA", "GATCGCGCCAA", "GATCGCGCCAA"]
//!     .iter()
//!     .map(|seq| Query::new(backbone.to_string(), seq.to_string(), 0, backbone.len()))
//!     .collect::<Vec<_>>();
//!
//! let consensus = sparc_consensus(backbone, &queries, &SparcConfig::default()).unwrap();
//! assert!(!consensus.seq.is_empty());
//! ```
//!
//! Inputs are validated before the algorithm runs; invalid input returns
//! [`SparcError`] instead of misbehaving.

mod best_path;
mod encoding;
mod error;
mod graph;
mod m5;
mod normalize;
mod pipeline;

pub use error::SparcError;
pub use m5::parse_m5;

/// 算法参数。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SparcConfig {
    /// debug=true 时向进程当前目录写调试文件
    /// （align_profile.txt / subgraph.dot / DEBUG.consensus.fasta 等）。
    pub debug: bool,
    /// k-mer 大小，范围 [1, 16]（k-mer 以 uint32 位编码，2 bit/碱基）。建议 [1, 2]。
    pub kmer: i32,
    /// 覆盖度阈值（CLI 的 c），建议 [1, 5]。
    pub coverage_threshold: i32,
    /// 打分方法：1 = 对数比例法，2 = 默认的线性减法。
    pub scoring_method: i32,
    /// debug 模式输出子图 dot 文件的区间起点（其余情况下不生效）。
    pub subgraph_begin: i32,
    /// debug 模式输出子图 dot 文件的区间终点（其余情况下不生效）。
    pub subgraph_end: i32,
    /// 保留字段，当前实现未使用。
    pub cns_start: i32,
    /// 保留字段，当前实现未使用。
    pub cns_end: i32,
    /// 保留字段：仅对原 CLI 读入的 m5 行生效，不影响通过 API 传入的 query。
    pub report_begin: i32,
    /// 保留字段：仅对原 CLI 读入的 m5 行生效，不影响通过 API 传入的 query。
    pub report_end: i32,
    /// 覆盖度滑动窗口半径（原 CLI 固定 200，绑定默认 2）。
    pub cov_radius: i32,
    /// 自适应阈值（CLI 的 t）。<0 关闭自适应（CLI 默认 -0.1），建议 [0.0, 0.3]。
    pub threshold: f64,
}

impl Default for SparcConfig {
    fn default() -> Self {
        Self {
            debug: false,
            kmer: 1,
            coverage_threshold: 2,
            scoring_method: 2,
            subgraph_begin: 0,
            subgraph_end: 0,
            cns_start: 0,
            cns_end: 0,
            report_begin: 0,
            report_end: 0,
            cov_radius: 2,
            threshold: 0.2,
        }
    }
}

/// consensus 结果。
///
/// `start`/`end` 是最优路径在 backbone 上的节点下标区间 `[start, end)`：
/// - `start == None`：路径头部不在 backbone 上（例如起点是一个插入分支节点）；
/// - `end == 0`：无可信路径（如无 query 输入），此时 `seq` 为原 backbone。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Consensus {
    pub seq: String,
    pub start: Option<u32>,
    pub end: u32,
}

/// 一条 read 相对 backbone 的比对，对应一行 blasr m5 记录。
///
/// `query_aligned_sequence` / `target_aligned_sequence` 是带 `-` gap 列的比对串，
/// 两者必须等长；`target_aligned_sequence` 中非 `-` 的碱基数必须等于
/// `target_end - target_start`。以上约束由 [`sparc_consensus`] 在调用时校验。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Query {
    pub(crate) query_aligned_sequence: String,
    pub(crate) target_aligned_sequence: String,
    pub(crate) is_reverse_strand: bool,
    pub(crate) query_start: usize,
    pub(crate) query_end: usize,
    pub(crate) target_start: usize,
    pub(crate) target_end: usize,
}

impl Query {
    /// 由 target/query 两条比对串构造（参数顺序：先 target 后 query）。
    /// `target_aligned` 为 backbone 方向的比对串，坐标为 backbone 上的
    /// `[target_start, target_end)`。默认按正链（tStrand = '+'）处理。
    pub fn new(
        target_aligned: String,
        query_aligned: String,
        target_start: usize,
        target_end: usize,
    ) -> Self {
        Self {
            query_aligned_sequence: query_aligned,
            target_aligned_sequence: target_aligned,
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
}

/// ConsensusNode.kmer 为 uint32_t，2 bit/碱基，最多编码 16 个碱基。
const MAX_KMER_LENGTH: i32 = 16;

/// 是否为标准碱基 ACGT（大小写不敏感），不含 gap。
fn is_standard_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

fn validate_inputs(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<(), SparcError> {
    if config.kmer < 1 || config.kmer > MAX_KMER_LENGTH {
        return Err(SparcError::InvalidKmer { kmer: config.kmer });
    }
    if config.cov_radius < 0 {
        return Err(SparcError::InvalidCovRadius {
            cov_radius: config.cov_radius,
        });
    }
    if config.threshold.is_nan() {
        return Err(SparcError::InvalidThreshold {
            threshold: config.threshold,
        });
    }
    if backbone.len() < config.kmer as usize {
        return Err(SparcError::BackboneTooShort {
            backbone_len: backbone.len(),
            kmer: config.kmer,
        });
    }
    for (position, &base) in backbone.as_bytes().iter().enumerate() {
        if !is_standard_base(base) {
            return Err(SparcError::InvalidBackboneBase {
                base: base as char,
                position,
            });
        }
    }
    for query in queries {
        validate_query_alignment(query, backbone.len())?;
    }
    Ok(())
}

fn validate_query_alignment(query: &Query, backbone_len: usize) -> Result<(), SparcError> {
    let query_bytes = query.query_aligned_sequence.as_bytes();
    let target_bytes = query.target_aligned_sequence.as_bytes();
    if query_bytes.is_empty() || target_bytes.is_empty() {
        return Err(SparcError::EmptyAlignment);
    }
    if query_bytes.len() != target_bytes.len() {
        return Err(SparcError::AlignedLengthMismatch {
            query_len: query_bytes.len(),
            target_len: target_bytes.len(),
        });
    }
    for (position, &base) in query_bytes.iter().enumerate() {
        if !is_standard_base(base) && base != b'-' {
            return Err(SparcError::InvalidAlignedBase {
                base: base as char,
                position,
            });
        }
    }
    let mut target_base_count = 0usize;
    for (position, &base) in target_bytes.iter().enumerate() {
        if !is_standard_base(base) && base != b'-' {
            return Err(SparcError::InvalidAlignedBase {
                base: base as char,
                position,
            });
        }
        if base != b'-' {
            target_base_count += 1;
        }
    }
    if query.target_start > query.target_end || query.target_end > backbone_len {
        return Err(SparcError::InvalidTargetSpan {
            target_start: query.target_start,
            target_end: query.target_end,
            backbone_len,
        });
    }
    if target_base_count != query.target_end - query.target_start {
        return Err(SparcError::SpanMismatch {
            span: query.target_end - query.target_start,
            target_bases: target_base_count,
        });
    }
    Ok(())
}

/// 由 backbone 和一组比对计算 consensus。
///
/// 返回 [`Consensus`]，其 `start`/`end` 为最优路径在 backbone 上的
/// `[start, end)` 节点下标区间（`end == 0` 表示回退为原 backbone，
/// `start == None` 表示路径头部不在 backbone 上）。
///
/// 所有输入在进入算法前完成校验；不满足约束时返回 [`SparcError`]。
pub fn sparc_consensus(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<Consensus, SparcError> {
    validate_inputs(backbone, queries, config)?;
    Ok(pipeline::compute_consensus(backbone, queries, config))
}

#[cfg(test)]
mod tests {

    use super::*;

    fn make_config(backbone_len: usize) -> SparcConfig {
        SparcConfig {
            debug: false,
            report_end: backbone_len as i32,
            subgraph_end: backbone_len as i32,
            cns_end: backbone_len as i32,
            ..SparcConfig::default()
        }
    }

    fn make_queries(backbone: &str) -> Vec<Query> {
        let q = |query_aligned_sequence: &str| Query {
            query_aligned_sequence: query_aligned_sequence.to_string(),
            target_aligned_sequence: backbone.to_string(),
            is_reverse_strand: false,
            query_start: 0,
            query_end: query_aligned_sequence.len(),
            target_start: 0,
            target_end: backbone.len(),
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
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
        assert_eq!(consensus.start, Some(0));
        assert!(consensus.end > 0);
    }

    /// 回归：空 queries 时必须走 fallback 路径，原样返回 backbone
    /// （seq == backbone、start == None、end == 0）。
    #[test]
    fn test_empty_queries_returns_backbone() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());

        let consensus = sparc_consensus(backbone, &[], &config).unwrap();
        assert_eq!(consensus.seq, backbone);
        assert_eq!(consensus.start, None);
        assert_eq!(consensus.end, 0);
    }

    /// 回归：重复调用 50 次结果稳定，不存在跨调用残留状态。
    #[test]
    fn test_sparc_consensus_repeated() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries = make_queries(backbone);

        for _ in 0..50 {
            let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
            assert!(!consensus.seq.is_empty());
            assert!(consensus.start.expect("start should be set") < consensus.end);
        }
    }

    /// 回归：3' 端插入使非 backbone 节点挂在最后一个 backbone 节点上，
    /// 该分支路径必须参与最优路径计算且不破坏结果。
    #[test]
    fn test_consensus_with_terminal_insertion() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let mut queries = make_queries(backbone);

        // target 比对串末尾为 '-'，即在 backbone 末端插入一个碱基
        let ins = |query_aligned_sequence: &str| Query {
            query_aligned_sequence: query_aligned_sequence.to_string(),
            target_aligned_sequence: format!("{backbone}-"),
            is_reverse_strand: false,
            query_start: 0,
            query_end: query_aligned_sequence.len(),
            target_start: 0,
            target_end: backbone.len(),
        };
        queries.push(ins("GATCGGGCTAAA"));
        queries.push(ins("GATCGGGCTAAG"));

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
        assert!(consensus.start.expect("start should be set") < consensus.end);
    }

    /// 回归：比对起始处即 mismatch 的 read 会在链首产生孤儿节点
    /// （不挂在任何 backbone 节点的右子图内），结果仍须非空。
    #[test]
    fn test_leading_mismatch_queries() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let queries: Vec<Query> = (0..1000)
            .map(|_| Query {
                query_aligned_sequence: "TATCGGGCTAA".to_string(),
                target_aligned_sequence: backbone.to_string(),
                is_reverse_strand: false,
                query_start: 0,
                query_end: 11,
                target_start: 0,
                target_end: backbone.len(),
            })
            .collect();

        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert!(!consensus.seq.is_empty());
    }

    /// 负链：m5 负链行的两条比对串按 read 方向给出，算法内部反向互补后
    /// 应与直接给正向比对得到完全一致的结果。
    #[test]
    fn test_reverse_strand_matches_forward() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());

        let forward = make_queries(backbone);
        let reversed: Vec<Query> = forward
            .iter()
            .map(|q| Query {
                query_aligned_sequence: reverse_complement(&q.query_aligned_sequence),
                target_aligned_sequence: reverse_complement(&q.target_aligned_sequence),
                is_reverse_strand: true,
                query_start: q.query_start,
                query_end: q.query_end,
                target_start: q.target_start,
                target_end: q.target_end,
            })
            .collect();

        let a = sparc_consensus(backbone, &forward, &config).unwrap();
        let b = sparc_consensus(backbone, &reversed, &config).unwrap();
        assert_eq!(a, b);
    }

    /// threshold 行为回归：自适应阈值必须真正影响结果。
    /// threshold=100 使所有边的 new_score 为负、从不更新，全部节点 score=0，
    /// 触发 fallback（原样返回 backbone）；threshold=-0.1 则正常产出 consensus。
    #[test]
    fn test_threshold_changes_behavior() {
        let backbone = "GATCGGGCTAA";
        let queries = make_queries(backbone);

        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = -0.1;
        let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
        assert_ne!(consensus.seq, backbone);
        assert!(consensus.end > 0);

        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = 100.0;
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
        let mut config = make_config(bb.len());
        config.kmer = 3;
        config.cov_radius = 500;
        let exact = Query {
            query_aligned_sequence: bb.to_string(),
            target_aligned_sequence: bb.to_string(),
            is_reverse_strand: false,
            query_start: 0,
            query_end: bb.len(),
            target_start: 0,
            target_end: bb.len(),
        };
        let mut queries = vec![exact.clone(), exact.clone(), exact];
        // 首列替换：产生孤儿链首；末端插入：产生末端分支
        queries.push(Query {
            query_aligned_sequence: format!("T{}", &bb[1..]),
            target_aligned_sequence: bb.to_string(),
            is_reverse_strand: false,
            query_start: 0,
            query_end: bb.len(),
            target_start: 0,
            target_end: bb.len(),
        });
        queries.push(Query {
            query_aligned_sequence: format!("{bb}A"),
            target_aligned_sequence: format!("{bb}-"),
            is_reverse_strand: false,
            query_start: 0,
            query_end: bb.len() + 1,
            target_start: 0,
            target_end: bb.len(),
        });

        let first = sparc_consensus(bb, &queries, &config).unwrap();
        assert!(!first.seq.is_empty());
        // 确定性：同输入重复调用结果一致
        for _ in 0..3 {
            let again = sparc_consensus(bb, &queries, &config).unwrap();
            assert_eq!(again, first);
        }
    }

    #[test]
    fn test_validation_errors() {
        let backbone = "GATCGGGCTAA";
        let config = make_config(backbone.len());
        let q =
            |query_aligned: &str, target_aligned: &str, span_start: usize, span_end: usize| Query {
                query_aligned_sequence: query_aligned.to_string(),
                target_aligned_sequence: target_aligned.to_string(),
                is_reverse_strand: false,
                query_start: 0,
                query_end: query_aligned.len(),
                target_start: span_start,
                target_end: span_end,
            };

        // kmer 越界（下界 / 上界）
        for kmer in [0, 17] {
            let mut cfg = make_config(backbone.len());
            cfg.kmer = kmer;
            assert_eq!(
                sparc_consensus(backbone, &[], &cfg),
                Err(SparcError::InvalidKmer { kmer })
            );
        }

        // backbone 长度不足
        let mut cfg = make_config(3);
        cfg.kmer = 5;
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

        // cov_radius 为负
        let mut cfg = make_config(backbone.len());
        cfg.cov_radius = -1;
        assert_eq!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidCovRadius { cov_radius: -1 })
        );

        // threshold 为 NaN（NaN != NaN，用 matches! 而非 assert_eq! 比较）
        let mut cfg = make_config(backbone.len());
        cfg.threshold = f64::NAN;
        assert!(matches!(
            sparc_consensus(backbone, &[], &cfg),
            Err(SparcError::InvalidThreshold { .. })
        ));

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
        assert_eq!(q.query_aligned_sequence, "GATCGCGCTAA");
        assert_eq!(q.target_aligned_sequence, "GATCGGGCTAA");
        assert_eq!(q.query_start, 0);
        assert_eq!(q.query_end, 12);
        assert_eq!(q.target_start, 0);
        assert_eq!(q.target_end, 11);
        assert!(!q.is_reverse_strand);

        let negative_row = "read_24/0_10 10 0 12 + ref_template 11 0 11 - 100 10 0 0 0 60 \
                            GATCGCGCTAA 10M GATCGGGCTAA";
        let q = Query::from_m5_row(negative_row).unwrap();
        assert!(q.is_reverse_strand);

        let bad = Query::from_m5_row("only three fields");
        assert!(matches!(bad, Err(SparcError::M5Parse(_))));

        let bad_strand = Query::from_m5_row(
            "r 10 0 12 + t 11 0 11 ? 100 10 0 0 0 60 GATCGCGCTAA 10M GATCGGGCTAA",
        );
        assert!(matches!(bad_strand, Err(SparcError::M5Parse(_))));
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
                .all(|q| q.target_start == 0 && !q.is_reverse_strand)
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
        let mut config = make_config(backbone.len());
        config.kmer = 1;
        config.threshold = -0.1;
        config.cov_radius = 200;

        let consensus = sparc_consensus(&backbone, &queries, &config).unwrap();
        let expected =
            std::fs::read_to_string(format!("{dir}/testdata2/out.consensus.fasta")).unwrap();
        let expected_seq: String = expected.lines().skip(1).collect();
        assert_eq!(consensus.seq, expected_seq);
    }
}
