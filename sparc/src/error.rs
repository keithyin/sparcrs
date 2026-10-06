//! [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。

/// [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum SparcError {
    /// kmer 必须 >= 1，且 <= 16（k-mer 以 uint32 位编码，2 bit/碱基）。
    #[error("invalid kmer {kmer}: expected 1..=16")]
    InvalidKmer { kmer: i32 },
    /// backbone 长度不足以切出 k-mer。
    #[error("backbone length {backbone_len} is shorter than kmer {kmer}")]
    BackboneTooShort { backbone_len: usize, kmer: i32 },
    /// backbone 含非法碱基（仅允许 ACGTacgt；含 N 时算法会静默出错，这里显式拒绝）。
    #[error("invalid base {base:?} in backbone at position {position}: expected ACGT")]
    InvalidBackboneBase { base: char, position: usize },
    /// cov_radius 为负会使滑动窗口下标越界。
    #[error("invalid cov_radius {cov_radius}: expected >= 0")]
    InvalidCovRadius { cov_radius: i32 },
    /// threshold 为 NaN。
    #[error("invalid threshold {threshold}: must not be NaN")]
    InvalidThreshold { threshold: f64 },
    /// 比对串为空。
    #[error("aligned sequences must not be empty")]
    EmptyAlignment,
    /// 两条比对串必须等长。
    #[error("aligned sequence length mismatch: query {query_len} != target {target_len}")]
    AlignedLengthMismatch { query_len: usize, target_len: usize },
    /// 比对串含非法字符（仅允许 ACGTacgt 与 gap '-'）。
    #[error("invalid base {base:?} in aligned sequence at position {position}: expected ACGT or '-'")]
    InvalidAlignedBase { base: char, position: usize },
    /// target 坐标超出 backbone 范围或 start > end。
    #[error("invalid target span [{target_start}, {target_end}): backbone length is {backbone_len}")]
    InvalidTargetSpan {
        target_start: usize,
        target_end: usize,
        backbone_len: usize,
    },
    /// target 比对串中非 gap 碱基数与声明的 span 不一致。
    #[error("target aligned sequence has {target_bases} non-gap bases, but span is {span}")]
    SpanMismatch { span: usize, target_bases: usize },
    /// m5 行解析失败。
    #[error("m5 parse error: {0}")]
    M5Parse(String),
}
