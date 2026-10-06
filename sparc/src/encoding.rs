//! 碱基的 2-bit 编码与 k-mer 位表示。
//!
//! 编码约定与原 C++ `str2bitsarr`/`bitsarr2str` 完全一致：
//! A/a→0、C/c→1、G/g→2、T/t→3，每个碱基 2 bit，
//! **第一个碱基位于最高位**（即 k-mer 值 = code_0 << 2(K-1) | ... | code_{K-1}）。
//! K ≤ 16 时 k-mer 值不超过 32 bit，可安全截断为 u32（对应 C++ 的
//! `current_node->kmer = (uint32_t)seq`）。

/// 单个碱基 → 2-bit 码。输入已由校验保证为 ACGTacgt。
pub(crate) fn base_to_two_bit_code(base: u8) -> u64 {
    debug_assert!(matches!(
        base,
        b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'
    ));
    match base {
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// 2-bit 码 → 碱基字符（大写）。
pub(crate) fn two_bit_code_to_base(code: u64) -> u8 {
    match code {
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'A',
    }
}

/// 将 2-bit 编码的 k-mer 解码为正向碱基序列（等价于 C++ `bitsarr2str`）。
pub(crate) fn decode_kmer_to_string(kmer: u32, kmer_size: usize) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(kmer_size);
    for offset in 0..kmer_size {
        let shift = 2 * (kmer_size - 1 - offset);
        decoded.push(two_bit_code_to_base(((kmer as u64) >> shift) & 3));
    }
    decoded
}

/// 原位反向互补。'-'（gap）与未知字符保持不变，等价于 C++ `reverse_complement_str`。
pub(crate) fn reverse_complement_in_place(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence.iter_mut() {
        *base = match *base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 与 C++ str2bitsarr 一致：首碱基在高位。
    #[test]
    fn test_encode_matches_str2bitsarr_semantics() {
        // "ACG" → (0 << 4) | (1 << 2) | 2 = 6
        let mut encoded = 0u64;
        for base in b"ACG" {
            encoded = (encoded << 2) | base_to_two_bit_code(*base);
        }
        assert_eq!(encoded, 6);
    }

    #[test]
    fn test_decode_kmer_roundtrip() {
        for sequence in ["A", "AC", "ACG", "ACGT", "ACGTACGTACGTACGT"] {
            let kmer_size = sequence.len();
            let mut encoded = 0u64;
            for base in sequence.bytes() {
                encoded = (encoded << 2) | base_to_two_bit_code(base);
            }
            let decoded = decode_kmer_to_string(encoded as u32, kmer_size);
            assert_eq!(decoded, sequence.as_bytes());
        }
    }

    #[test]
    fn test_reverse_complement_in_place() {
        // 反向互补：'-'（gap）保持不变
        let mut sequence = b"ACGT-N".to_vec();
        reverse_complement_in_place(&mut sequence);
        // 反转得 "N-TGCA"，互补得 "N-ACGT"
        assert_eq!(sequence, b"N-ACGT");
    }

    #[test]
    fn test_reverse_complement_is_involution() {
        let original = b"ACGTacgt-".to_vec();
        let mut sequence = original.clone();
        reverse_complement_in_place(&mut sequence);
        reverse_complement_in_place(&mut sequence);
        // 小写会被大写化，其余应复原
        assert_eq!(sequence, b"ACGTACGT-".to_vec());
        let _ = original;
    }
}
