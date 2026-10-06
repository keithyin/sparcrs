//! [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。

/// [`sparc_consensus`](crate::sparc_consensus) 与 m5 解析的错误类型。
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum SparcError {
    /// kmer 必须 >= 1，且 <= 16（k-mer 以 uint32 位编码，2 bit/碱基）。
    #[error("invalid kmer {kmer}: expected 1..=16")]
    InvalidKmer {
        /// 传入的 kmer 值。
        kmer: u8,
    },
    /// backbone 长度不足以切出 k-mer。
    #[error("backbone length {backbone_len} is shorter than kmer {kmer}")]
    BackboneTooShort {
        /// backbone 的字节数。
        backbone_len: usize,
        /// 传入的 kmer 值。
        kmer: u8,
    },
    /// backbone 含非法碱基（仅允许 ACGTacgt；含 N 时算法会静默出错，这里显式拒绝）。
    #[error("invalid base {base:?} in backbone at position {position}: expected ACGT")]
    InvalidBackboneBase {
        /// 非法碱基字符。
        base: char,
        /// 在 backbone 中的下标。
        position: usize,
    },
    /// threshold 为 NaN 或 +∞（NaN 会使所有比较恒假，+∞ 使自适应阈值失效）。
    #[error("invalid threshold {threshold}: must not be NaN or +infinity")]
    InvalidThreshold {
        /// 传入的 threshold 值。
        threshold: f64,
    },
    /// 比对串为空。
    #[error("aligned sequences must not be empty")]
    EmptyAlignment,
    /// 两条比对串必须等长。
    #[error("aligned sequence length mismatch: query {query_len} != target {target_len}")]
    AlignedLengthMismatch {
        /// query 比对串长度。
        query_len: usize,
        /// target 比对串长度。
        target_len: usize,
    },
    /// 比对串含非法字符（仅允许 ACGTacgt 与 gap '-'）。
    #[error(
        "invalid base {base:?} in aligned sequence at position {position}: expected ACGT or '-'"
    )]
    InvalidAlignedBase {
        /// 非法字符。
        base: char,
        /// 在比对串中的下标。
        position: usize,
    },
    /// target 坐标超出 backbone 范围或 start > end。
    #[error(
        "invalid target span [{target_start}, {target_end}): backbone length is {backbone_len}"
    )]
    InvalidTargetSpan {
        /// 声明的 span 起点。
        target_start: usize,
        /// 声明的 span 终点（不含）。
        target_end: usize,
        /// backbone 的字节数。
        backbone_len: usize,
    },
    /// target 比对串中非 gap 碱基数与声明的 span 不一致。
    #[error("target aligned sequence has {target_bases} non-gap bases, but span is {span}")]
    SpanMismatch {
        /// 声明的 span 长度。
        span: usize,
        /// 比对串中非 gap 的碱基数。
        target_bases: usize,
    },
    /// m5 行解析失败。`line` 为 1 起始的物理行号；单行解析（如
    /// [`Query::from_m5_row`](crate::Query::from_m5_row)）时为 1。
    #[error("m5 parse error at line {line}: {message}")]
    M5Parse {
        /// 1 起始的物理行号。
        line: usize,
        /// 具体解析失败原因。
        message: String,
    },
    /// debug 输出文件创建失败（仅 `SparcConfig::debug_output` 开启时可能发生）。
    /// 消息为 `<文件名>: <io 错误>`。
    #[error("failed to create debug output file: {0}")]
    DebugWrite(String),
}
