//! 新旧 Sparc 实现的差分对拍测试。
//!
//! 用确定性 RNG（splitmix64，固定种子，无外部依赖）生成随机比对场景，
//! 分别喂给纯 Rust 实现与 C++ binding（../Sparc），断言
//! `Consensus { seq, start, end }` 完全一致。
//!
//! 生成器刻意覆盖的 parity 关键路径：
//! - 插入/删除/替换混合的比对（触发归一化的 mismatch 拆分与 gap 交换）；
//! - 首列 mismatch（孤儿链首）；
//! - 短 backbone（node_vec 只有 1 个节点）；
//! - cov_radius 大于节点数（滑窗 clamp）、radius = 0（前缀最大值 quirk）；
//! - 负链比对（反向互补路径）；
//! - scoring_method 1（C++ 既有 bug 行为）/ 2 / 3；
//! - threshold 关闭自适应与极端大值（fallback 路径）。
//!
//! 已知边界：生成器不直接产生「双 gap 列」（同一列两条比对串都是 '-'），
//! 但 C++ 归一化的 gap 交换会在部分合法输入上制造出双 gap 列，进而触发
//! 其 match 区锚点越界的既有 UB（GraphConstruction.cpp:905
//! `node_vec[MatchPosition]` 越界，ASAN 实证），C++ 进程直接 SIGSEGV。
//! 因此对拍通过独立子进程 oracle 进行（见 `src/bin/legacy-oracle.rs`），
//! oracle 崩溃被记录为 C++ 既有 bug，不计入对拍失败。

use sparc::{Query, SparcConfig};
use sparc_legacy as legacy;

/// splitmix64 确定性随机数发生器
pub struct Random {
    state: u64,
}

impl Random {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// [0, bound) 内的随机下标；bound 必须非零
    pub fn below(&mut self, bound: usize) -> usize {
        (self.next_u64() % bound as u64) as usize
    }

    pub fn chance(&mut self, percent: u64) -> bool {
        self.next_u64() % 100 < percent
    }

    pub fn pick<T: Copy, const N: usize>(&mut self, choices: &[T; N]) -> T {
        choices[self.below(choices.len())]
    }
}

pub const BASES: &[u8] = b"ACGT";

pub fn random_base(rng: &mut Random) -> u8 {
    BASES[rng.below(4)]
}

pub fn reverse_complement(sequence: &str) -> String {
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

/// 一条随机比对：target 串为 backbone 切片原样（非 gap 数 == span），
/// query 串由 target 经替换/插入/删除构造。
#[derive(Clone)]
pub struct GeneratedAlignment {
    pub target_aligned: String,
    pub query_aligned: String,
    pub target_start: usize,
    pub target_end: usize,
    pub reverse_strand: bool,
}

pub fn generate_alignment(rng: &mut Random, backbone: &[u8]) -> GeneratedAlignment {
    let span_length = 1 + rng.below(backbone.len());
    generate_alignment_with_span(rng, backbone, span_length)
}

/// 指定 span 长度的比对生成（基准测试用，控制覆盖度）。
pub fn generate_alignment_with_span(
    rng: &mut Random,
    backbone: &[u8],
    span_length: usize,
) -> GeneratedAlignment {
    let span_length = span_length.min(backbone.len());
    let target_start = rng.below(backbone.len() - span_length + 1);
    let target_end = target_start + span_length;
    let target_sequence = &backbone[target_start..target_end];

    let mut target_columns: Vec<u8> = Vec::with_capacity(span_length + 8);
    let mut query_columns: Vec<u8> = Vec::with_capacity(span_length + 8);

    for (offset, &reference_base) in target_sequence.iter().enumerate() {
        // 首列强制 mismatch 的概率（孤儿链首场景）
        let force_mismatch = offset == 0 && rng.chance(15);
        match rng.below(100) {
            0..=54 => {
                // 匹配列
                let query_base = if force_mismatch {
                    BASES
                        .iter()
                        .find(|&&b| b != reference_base)
                        .copied()
                        .unwrap()
                } else {
                    reference_base
                };
                query_columns.push(query_base);
                target_columns.push(reference_base);
            }
            55..=69 => {
                // 替换列
                let mutated = BASES
                    .iter()
                    .filter(|&&b| b != reference_base)
                    .nth(rng.below(3))
                    .copied()
                    .unwrap();
                query_columns.push(mutated);
                target_columns.push(reference_base);
            }
            70..=79 => {
                // 删除列：query 侧 gap（保持 query 非空即可）
                query_columns.push(b'-');
                target_columns.push(reference_base);
            }
            _ => {
                // 插入列：query 多出 1~2 个碱基（target 侧 gap）
                let insertions = 1 + rng.below(2);
                for _ in 0..insertions {
                    query_columns.push(random_base(rng));
                    target_columns.push(b'-');
                }
                query_columns.push(reference_base);
                target_columns.push(reference_base);
            }
        }
    }

    let reverse_strand = rng.chance(20);
    let target_aligned = String::from_utf8(target_columns).unwrap();
    let query_aligned = String::from_utf8(query_columns).unwrap();
    let (target_aligned, query_aligned) = if reverse_strand {
        (
            reverse_complement(&target_aligned),
            reverse_complement(&query_aligned),
        )
    } else {
        (target_aligned, query_aligned)
    };

    GeneratedAlignment {
        target_aligned,
        query_aligned,
        target_start,
        target_end,
        reverse_strand,
    }
}

pub fn generate_config(rng: &mut Random, backbone_len: usize) -> SparcConfig {
    let kmer = rng.pick(&[1, 1, 2, 2, 3, 4]);
    let cov_radius = rng.pick(&[0, 1, 2, 3, 10, 200, 500, 1000]);
    let threshold = rng.pick(&[-1.0_f64, -0.1, 0.0, 0.1, 0.2, 0.3, 100.0]);
    let scoring_method = rng.pick(&[2, 2, 2, 1, 3]);
    let coverage_threshold = rng.pick(&[1, 2, 5]);
    SparcConfig {
        debug: false,
        kmer,
        coverage_threshold,
        scoring_method,
        subgraph_begin: 0,
        subgraph_end: backbone_len as i32,
        cns_start: 0,
        cns_end: backbone_len as i32,
        report_begin: 0,
        report_end: backbone_len as i32,
        cov_radius,
        threshold,
    }
}

/// 新 crate 的 config → 旧 crate 的 config（字段逐一对应）。
pub fn to_legacy_config(config: &SparcConfig) -> legacy::SparcConfig {
    legacy::SparcConfig {
        debug: config.debug,
        kmer: config.kmer,
        coverage_threshold: config.coverage_threshold,
        scoring_method: config.scoring_method,
        subgraph_begin: config.subgraph_begin,
        subgraph_end: config.subgraph_end,
        cns_start: config.cns_start,
        cns_end: config.cns_end,
        report_begin: config.report_begin,
        report_end: config.report_end,
        cov_radius: config.cov_radius,
        threshold: config.threshold,
    }
}

/// 生成一个随机场景（backbone + reads + config），并同时构造两套 Query。
pub struct Scenario {
    pub backbone: String,
    pub config: SparcConfig,
    pub legacy_config: legacy::SparcConfig,
    pub queries_new: Vec<Query>,
    pub queries_legacy: Vec<legacy::Query>,
    /// 保留原始比对数据（oracle 传输与调试用）
    pub raw_alignments: Vec<GeneratedAlignment>,
}

pub fn generate_scenario(rng: &mut Random) -> Scenario {
    let backbone_length = match rng.below(10) {
        0..=6 => 10 + rng.below(190),
        _ => 200 + rng.below(1800),
    };
    let mut backbone_bytes = Vec::with_capacity(backbone_length);
    for _ in 0..backbone_length {
        backbone_bytes.push(random_base(rng));
    }
    let backbone = String::from_utf8(backbone_bytes).unwrap();

    let read_count = match rng.below(10) {
        0 => 0, // 空 queries：fallback 路径
        1..=7 => 1 + rng.below(40),
        _ => 40 + rng.below(160),
    };

    let mut queries_new = Vec::with_capacity(read_count);
    let mut queries_legacy = Vec::with_capacity(read_count);
    let mut raw_alignments = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let alignment = generate_alignment(rng, backbone.as_bytes());
        raw_alignments.push(GeneratedAlignment {
            target_aligned: alignment.target_aligned.clone(),
            query_aligned: alignment.query_aligned.clone(),
            target_start: alignment.target_start,
            target_end: alignment.target_end,
            reverse_strand: alignment.reverse_strand,
        });
        let new_query = Query::new(
            alignment.target_aligned.clone(),
            alignment.query_aligned.clone(),
            alignment.target_start,
            alignment.target_end,
        );
        let new_query = if alignment.reverse_strand {
            new_query.reverse_strand()
        } else {
            new_query
        };
        let legacy_query = legacy::Query::new(
            alignment.target_aligned,
            alignment.query_aligned,
            alignment.target_start,
            alignment.target_end,
        );
        let legacy_query = if alignment.reverse_strand {
            legacy_query.reverse_strand()
        } else {
            legacy_query
        };
        queries_new.push(new_query);
        queries_legacy.push(legacy_query);
    }

    // 刻意场景：偶尔整体重复一份 reads（重复 kmer / 边覆盖度叠加）
    if rng.chance(10) && !queries_new.is_empty() {
        raw_alignments.extend(raw_alignments.clone());
        queries_new.extend(queries_new.clone());
        queries_legacy.extend(queries_legacy.clone());
    }

    let config = generate_config(rng, backbone.len());
    let legacy_config = to_legacy_config(&config);
    Scenario {
        backbone,
        config,
        legacy_config,
        queries_new,
        queries_legacy,
        raw_alignments,
    }
}

pub fn serialize_scenario(scenario: &Scenario) -> String {
    let mut fields: Vec<String> = vec![
        scenario.backbone.clone(),
        scenario.config.kmer.to_string(),
        scenario.config.coverage_threshold.to_string(),
        scenario.config.scoring_method.to_string(),
        scenario.config.subgraph_begin.to_string(),
        scenario.config.subgraph_end.to_string(),
        scenario.config.cns_start.to_string(),
        scenario.config.cns_end.to_string(),
        scenario.config.report_begin.to_string(),
        scenario.config.report_end.to_string(),
        scenario.config.cov_radius.to_string(),
        scenario.config.threshold.to_string(),
        scenario.raw_alignments.len().to_string(),
    ];
    for alignment in &scenario.raw_alignments {
        fields.push(alignment.target_start.to_string());
        fields.push(alignment.target_end.to_string());
        fields.push(if alignment.reverse_strand {
            "1".into()
        } else {
            "0".into()
        });
        fields.push(alignment.query_aligned.clone());
        fields.push(alignment.target_aligned.clone());
    }
    fields.join("\t")
}

/// oracle 的结果：成功时为 `(seq, start, end)`，失败时为错误消息。
pub type OracleResult = Result<(String, Option<u32>, u32), String>;

/// 解析 oracle 的输出行：`OK \t start \t end \t seq` 或 `ERR \t message`。
pub fn parse_oracle_line(line: &str) -> Result<OracleResult, String> {
    let fields: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
    match fields[0] {
        "OK" => {
            let start = if fields[1] == "-1" {
                None
            } else {
                Some(fields[1].parse::<u32>().expect("oracle start 解析失败"))
            };
            let end = fields[2].parse::<u32>().expect("oracle end 解析失败");
            Ok(Ok((fields[3].to_string(), start, end)))
        }
        "ERR" => Ok(Err(fields[1].to_string())),
        other => Err(format!("oracle 输出无法解析: {other:?}")),
    }
}
