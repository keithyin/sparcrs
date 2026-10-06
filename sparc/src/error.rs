//! [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。

use std::fmt;

/// [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。
#[derive(Clone, Debug, PartialEq)]
pub enum SparcError {
    /// kmer 必须 >= 1，且 <= 16（k-mer 以 uint32 位编码，2 bit/碱基）。
    InvalidKmer { kmer: i32 },
    /// backbone 长度不足以切出 k-mer。
    BackboneTooShort { backbone_len: usize, kmer: i32 },
    /// backbone 含非法碱基（仅允许 ACGTacgt；含 N 时算法会静默出错，这里显式拒绝）。
    InvalidBackboneBase { base: char, position: usize },
    /// cov_radius 为负会使滑动窗口下标越界。
    InvalidCovRadius { cov_radius: i32 },
    /// threshold 为 NaN。
    InvalidThreshold { threshold: f64 },
    /// 比对串为空。
    EmptyAlignment,
    /// 两条比对串必须等长。
    AlignedLengthMismatch { query_len: usize, target_len: usize },
    /// 比对串含非法字符（仅允许 ACGTacgt 与 gap '-'）。
    InvalidAlignedBase { base: char, position: usize },
    /// target 坐标超出 backbone 范围或 start > end。
    InvalidTargetSpan {
        target_start: usize,
        target_end: usize,
        backbone_len: usize,
    },
    /// target 比对串中非 gap 碱基数与声明的 span 不一致。
    SpanMismatch { span: usize, target_bases: usize },
    /// m5 行解析失败。
    M5Parse(String),
}

impl fmt::Display for SparcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SparcError::InvalidKmer { kmer } => {
                write!(f, "invalid kmer {kmer}: expected 1..=16")
            }
            SparcError::BackboneTooShort { backbone_len, kmer } => {
                write!(
                    f,
                    "backbone length {backbone_len} is shorter than kmer {kmer}"
                )
            }
            SparcError::InvalidBackboneBase { base, position } => {
                write!(
                    f,
                    "invalid base {base:?} in backbone at position {position}: expected ACGT"
                )
            }
            SparcError::InvalidCovRadius { cov_radius } => {
                write!(f, "invalid cov_radius {cov_radius}: expected >= 0")
            }
            SparcError::InvalidThreshold { threshold } => {
                write!(f, "invalid threshold {threshold}: must not be NaN")
            }
            SparcError::EmptyAlignment => write!(f, "aligned sequences must not be empty"),
            SparcError::AlignedLengthMismatch {
                query_len,
                target_len,
            } => {
                write!(
                    f,
                    "aligned sequence length mismatch: query {query_len} != target {target_len}"
                )
            }
            SparcError::InvalidAlignedBase { base, position } => {
                write!(
                    f,
                    "invalid base {base:?} in aligned sequence at position {position}: expected ACGT or '-'"
                )
            }
            SparcError::InvalidTargetSpan {
                target_start,
                target_end,
                backbone_len,
            } => {
                write!(
                    f,
                    "invalid target span [{target_start}, {target_end}): backbone length is {backbone_len}"
                )
            }
            SparcError::SpanMismatch { span, target_bases } => {
                write!(
                    f,
                    "target aligned sequence has {target_bases} non-gap bases, but span is {span}"
                )
            }
            SparcError::M5Parse(msg) => write!(f, "m5 parse error: {msg}"),
        }
    }
}

impl std::error::Error for SparcError {}
