//! 输入校验：在进入算法前完成，不满足约束返回 [`SparcError`]。
//!
//! C++ 绑定层没有这些防护（非法输入直接触发未定义行为），这里集中补齐；
//! 校验只拒绝算法无法正确定义的输入，不改变合法输入的计算语义。

use crate::Query;
use crate::config::SparcConfig;
use crate::error::SparcError;

/// ConsensusNode.kmer 为 uint32_t，2 bit/碱基，最多编码 16 个碱基。
pub(crate) const MAX_KMER_LENGTH: u8 = 16;

/// 是否为标准碱基 ACGT（大小写不敏感），不含 gap。
pub(crate) fn is_standard_base(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

pub(crate) fn validate_inputs(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Result<(), SparcError> {
    if config.kmer == 0 || config.kmer > MAX_KMER_LENGTH {
        return Err(SparcError::InvalidKmer { kmer: config.kmer });
    }
    if config.threshold.is_nan() || config.threshold == f64::INFINITY {
        return Err(SparcError::InvalidThreshold {
            threshold: config.threshold,
        });
    }
    if backbone.len() < usize::from(config.kmer) {
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
