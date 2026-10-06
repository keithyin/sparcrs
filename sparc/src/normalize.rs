//! 比对归一化（对应 C++ `SparcNormalizeAlignment`）。
//!
//! 将 mismatch 列拆成一个"删除列 + 插入列"对（gap 靠左/右的排布与 C++ 一致），
//! 并在扫描过程中原位改写输入串（前移可配对的碱基、把原来的位置置为 gap）。
//! 槽位推进逻辑与 C++ 逐行对应，包括 mismatch 分支里反直觉的
//! `n_char--` 回退写法。

/// 归一化两条等长的比对串（就地替换为新的两条串）。
///
/// 输入串仅含 ACGTacgt 与 '-'（已由输入校验保证），长度相等。
pub(crate) fn normalize_alignment(query_aligned: &mut Vec<u8>, target_aligned: &mut Vec<u8>) {
    let sequence_length = query_aligned.len();
    // mismatch 列 1 → 2，最坏情况下新串长度为原来的 2 倍
    let mut query_new = vec![0u8; 2 * sequence_length];
    let mut target_new = vec![0u8; 2 * sequence_length];
    let mut chars_written = 0usize;

    for column in 0..sequence_length {
        let query_base = query_aligned[column];
        let target_base = target_aligned[column];

        if query_base == target_base {
            // 匹配列（包括双 gap 列）原样保留
            query_new[chars_written] = query_base;
            target_new[chars_written] = target_base;
            chars_written += 1;
        } else if query_base != b'-' && target_base != b'-' {
            // mismatch 列：target 碱基先占住一列
            target_new[chars_written] = target_base;
            chars_written += 1;

            // 向右找 target 上第一个非 gap 碱基，若与 query 碱基相同，
            // 则把它挪走（置 gap），并用 query 碱基自身配对
            let mut replaced = false;
            if let Some(look_ahead) = position_of_first_non_gap(target_aligned, column) {
                if target_aligned[look_ahead] == query_base {
                    replaced = true;
                    target_aligned[look_ahead] = b'-';
                    target_new[chars_written] = query_base;
                }
            }
            if !replaced {
                target_new[chars_written] = b'-';
            }
            // C++ 在此回退 n_char：上一格写入的 q[i]/'-' 将与 q_new 的 '-' 同格
            chars_written -= 1;
            query_new[chars_written] = b'-';
            chars_written += 1;
            query_new[chars_written] = query_base;
            chars_written += 1;
        } else {
            // 单侧 gap 列：向右找同串第一个非 gap 碱基，
            // 若与对面碱基相同则交换（把 gap 右移一格）
            if query_base == b'-' {
                shift_gap_right(query_aligned, column, target_base);
            } else {
                shift_gap_right(target_aligned, column, query_base);
            }
            // 写入（可能已被交换过的）当前列
            query_new[chars_written] = query_aligned[column];
            target_new[chars_written] = target_aligned[column];
            chars_written += 1;
        }
    }

    query_new.truncate(chars_written);
    target_new.truncate(chars_written);
    *query_aligned = query_new;
    *target_aligned = target_new;
}

/// `sequence[range_start + 1..]` 中第一个非 gap 碱基的下标。
fn position_of_first_non_gap(sequence: &[u8], range_start: usize) -> Option<usize> {
    sequence[range_start + 1..]
        .iter()
        .position(|&base| base != b'-')
        .map(|offset| range_start + 1 + offset)
}

/// 把 `sequence[column]` 的 gap 与右侧第一个可配对碱基交换
/// （等价于 C++ 中"向右找同串首个非 gap 碱基，配对则交换"的循环：
/// 首个非 gap 碱基不配对时同样终止，不交换）。
fn shift_gap_right(sequence: &mut [u8], column: usize, partner_base: u8) {
    if let Some(look_ahead) = position_of_first_non_gap(sequence, column) {
        if sequence[look_ahead] == partner_base {
            sequence[column] = sequence[look_ahead];
            sequence[look_ahead] = b'-';
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn normalize(query: &str, target: &str) -> (String, String) {
        let mut q = query.as_bytes().to_vec();
        let mut t = target.as_bytes().to_vec();
        normalize_alignment(&mut q, &mut t);
        (String::from_utf8(q).unwrap(), String::from_utf8(t).unwrap())
    }

    #[test]
    fn test_identical_sequences_unchanged() {
        assert_eq!(normalize("ACGT", "ACGT"), ("ACGT".into(), "ACGT".into()));
    }

    /// mismatch 列拆成两列：target 碱基占一列（配 query 的 gap），
    /// query 碱基占一列（与 target 同碱基自配对，或配 gap）。
    #[test]
    fn test_mismatch_split_with_lookahead_match() {
        // q = AGT, t = ACT：列1 G vs C mismatch，
        // 向右找 t 首个非 gap：T ≠ G → 不替换 → q="A-GT", t="AC-T"
        assert_eq!(normalize("AGT", "ACT"), ("A-GT".into(), "AC-T".into()));
    }

    /// mismatch 列且 target 后续存在相同碱基：该碱基被置 gap，query 碱基自配对。
    #[test]
    fn test_mismatch_split_with_replacement() {
        // q = AGG, t = ACG：列1 G vs C mismatch，
        // t[2]='G' == q[1]='G' → t[2] 置 '-'，query G 与自身配对
        // → q="A-GG", t="ACG-"
        assert_eq!(normalize("AGG", "ACG"), ("A-GG".into(), "ACG-".into()));
    }

    /// 单侧 gap 列：向右找同串首个非 gap 碱基，配对则交换（gap 右移一格）。
    #[test]
    fn test_query_gap_swap() {
        // q = A-T, t = AT-：列1 q gap，q[2]='T' == t[1]='T' → 交换后
        // 列1 变成 (T,T) 匹配列，q[2] 变 '-'；列2 成双 gap 原样保留
        // → q="AT-", t="AT-"
        assert_eq!(normalize("A-T", "AT-"), ("AT-".into(), "AT-".into()));
    }

    /// 单侧 gap 列但后续碱基不配对：保持原样。
    #[test]
    fn test_query_gap_no_swap() {
        // q = A-G, t = AT-：列1 q gap，q[2]='G' ≠ t[1]='T' → 不交换
        // 列2 双 gap 原样 → 不变
        assert_eq!(normalize("A-G", "AT-"), ("A-G".into(), "AT-".into()));
    }

    /// target 侧 gap 的对称交换。
    #[test]
    fn test_target_gap_swap() {
        // q = AT-, t = A-T：列1 t gap，t[2]='T' == q[1]='T' → 交换后
        // 列1 变 (T,T)，t[2] 变 '-'；列2 双 gap 原样 → q="AT-", t="AT-"
        assert_eq!(normalize("AT-", "A-T"), ("AT-".into(), "AT-".into()));
    }

    #[test]
    fn test_mismatch_expands_length_by_one() {
        // 全 mismatch：每列 1 → 2，长度翻倍
        let (q, t) = normalize("AC", "TG");
        // 列0: T 占位、A 配 gap；列1: G 占位、C 配 gap（无后续可配碱基）
        assert_eq!(q.len(), 4);
        assert_eq!(t.len(), 4);
        assert_eq!(
            (q.as_bytes(), t.as_bytes()),
            (b"-A-C".as_slice(), b"T-G-".as_slice())
        );
    }
}
