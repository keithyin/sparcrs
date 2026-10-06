//! 比对归一化（对应 C++ `SparcNormalizeAlignment`）。
//!
//! 将 mismatch 列拆成一个"删除列 + 插入列"对（gap 靠左/右的排布与 C++ 一致），
//! 并在扫描过程中原位改写输入串（前移可配对的碱基、把原来的位置置为 gap）。
//! 槽位推进逻辑与 C++ 逐行对应，包括 mismatch 分支里反直觉的
//! `n_char--` 回退写法。

/// 归一化两条等长的比对串（就地替换为归一化后的新串），使下游图构建
/// （`graph.rs` 的 `add_path_to_backbone`）可以逐列无歧义地消费。
///
/// 对应 C++ `SparcNormalizeAlignment`：槽位推进逻辑（包括反直觉的游标回退）
/// 与其逐行对应，输出逐字节一致，由 parity_tests 差分测试锁定。
///
/// # 做了什么
///
/// `query_aligned`（read 侧）与 `target_aligned`（backbone 侧）只含
/// `ACGTacgt` 与 `'-'` 且等长。逐列扫描，按列的形态分三种处理：
///
/// 1. **相等列**（两侧相同，双 gap 列也算相等）：原样复制。
/// 2. **mismatch 列**（两侧均非 gap 但不等）：拆成两列。先写**删除列**
///    「target 碱基 + query gap」，再写**插入列**「query 碱基 + target gap」；
///    写插入列前向右看 target 串，若第一个非 gap 碱基（跳过连续 gap）恰好
///    等于 query 碱基，则把它在**输入串**中置为 gap（预支过来）与 query
///    碱基配对，插入列成为「query 碱基 + 同一碱基」的相等列。
/// 3. **单侧 gap 列**：在 gap 所在串内向右找第一个非 gap 碱基（跳过连续
///    gap），与对面串当前列的碱基相同则交换——gap 右移一格、当前列变成
///    相等列；不配对或已到串尾则原样写入。
///
/// 第 2、3 种处理都会**就地改写输入串**（预支/换走的碱基变成 gap），后续列
/// 的判定能看到这些改动；每列只推一步，跨列累积等效于把 gap 一路推过整段
/// 可配对的碱基。
///
/// 结果写进两条新缓冲区（按 2n 预留：每个 mismatch 列 1 → 2，最坏翻倍），
/// 截断到实际列数后替换原串。输出保证：
///
/// * 不再存在「两侧均非 gap 且互不相等」的列——每列要么相等、要么至少
///   一侧是 gap；
/// * 忽略 gap 后，两条串各自的碱基序列与输入相同（只移动了 gap 位置和
///   配对关系，不增删碱基）；
/// * 长度可能增长（见上），调用方需以替换后的串为准。
///
/// # 为什么这么做
///
/// 下游 `add_path_to_backbone` 对这对串的用法完全是列级的：
///
/// * query 侧为 gap 的列不产生节点，只推进 backbone 位置（target 非 gap 时）；
/// * 从某列起连续 K 个列两两相等（双 gap 列算相等）才算 match 区，把
///   query 侧 k-mer 锚定到 backbone 节点并 +1 覆盖度；
/// * 其余窗口属于分支区，从该列起收集 query 侧前 K 个非 gap 碱基编码为
///   分支 k-mer（target 侧碱基永远不会进入 k-mer）。
///
/// 在这套规则下，mismatch 列是歧义列：一列同时压着一个 read 没有的
/// backbone 碱基（应按"跳过、不计覆盖度"处理）和一个 backbone 没有的
/// read 碱基（应按"插入、进分支"处理）。拆成删除列 + 插入列后，每列只剩
/// 一种角色，列级规则对每列都有唯一解释。同理，比对器放置 gap 的位置有
/// 任意性（序列等价、排布不同），而 match 窗口按列判定——交换/预支把序列
/// 上成立、排布上被掩盖的配对找回来；否则与 backbone 完全一致的 read 也拿
/// 不到覆盖度，反而多出分支节点（见示例 3、4）。
///
/// # 示例
///
/// ```text
/// (1) mismatch 拆分（右侧没有可预支的相同碱基）
///
///     输入  q = A G T          输出  q = A - G T
///           t = A C T                t = A C - T
///               ^^^^^^^                    ^^^^  ^^^^
///               G/C 错配列                 删除列  插入列
///
///     错配被改写成：read 在该位置没有碱基对应 backbone 的 C（删除），
///     read 的 G 没有 backbone 碱基与之对应（插入）。
///
/// (2) mismatch 拆分 + 预支：t 右侧第一个非 gap 碱基 G 恰好等于 query 的 G
///
///     输入  q = A G G          输出  q = A - G G
///           t = A C G                t = A C G -
///               ^^^^^^^                    ^^^^  ^^^^  ^^^^
///               G/C 错配列                 删除列  配对列  插入列
///
///     t 的 G 被预支（原位置置 gap），read 的第一个 G 与它配对成 (G,G)
///     相等列；read 的第二个 G 落在插入列。若不预支，则是 read 的第一个
///     G 配 gap、第二个 G 配 backbone 的 G——碱基内容相同、配对关系不同，
///     预支的排布是 C++ 选定的规范形。
///
/// (3) 单侧 gap 交换：两条序列完全相同，仅 gap 排布不同
///
///     输入  q = A - T          输出  q = A T -
///           t = A T -                t = A T -
///
///     q 的 gap 与右侧的 T 交换：列 1 变成 (T,T) 相等列，原 gap 与 t 的
///     gap 拼成双 gap 列。K=2 时输入排布没有任何全等窗口，输出排布的
///     第一个窗口即全等——被排布掩盖的 match 被找回。
///
/// (4) 预支 + 交换联动（K=2，backbone = ACGG）
///
///     输入  q = A G G T        输出  q = A - G G T
///           t = A C G G              t = A C G G -
///
///     列 1 的 G/C 错配预支了 t 的第一个 G；列 2 变成单侧 gap 列后又把
///     t 的第二个 G 换上来，于是输出的列 2、3 都是 (G,G)，K=2 的全等
///     窗口把 read 的 "GG" 锚定到 backbone 的 "GG" 节点。不归一化时三个
///     窗口都含不等列，整个 read 只会生成 3 个分支节点、backbone 零覆盖。
/// ```
///
/// # 性能
///
/// 每列向右最多扫描到第一个非 gap 碱基；常规比对中 gap 很短，近似线性。
/// 最坏情况（长 gap 串配同碱基段，扫描步长逐列增长）为 O(n²)，与 C++ 相同。
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
            // 相等列（含双 gap 列）原样保留
            query_new[chars_written] = query_base;
            target_new[chars_written] = target_base;
            chars_written += 1;
        } else if query_base != b'-' && target_base != b'-' {
            // mismatch 列拆两列。第一列（删除列）：target 碱基占位、query 侧留 gap
            target_new[chars_written] = target_base;
            chars_written += 1;

            // 第二列（插入列）的 target 侧：向右找 target 上第一个非 gap 碱基，
            // 若与 query 碱基相同则预支——在输入串中置为 gap，与 query 碱基
            // 配对成相等列；否则插入列配 gap
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
            // C++ 在此回退 n_char：游标退回删除列补写 query 侧的 gap，再前进
            // 一格写插入列的 query 碱基——一个游标交错写完两列的四格
            chars_written -= 1;
            query_new[chars_written] = b'-';
            chars_written += 1;
            query_new[chars_written] = query_base;
            chars_written += 1;
        } else {
            // 单侧 gap 列：向右找 gap 所在串的第一个非 gap 碱基，与对面碱基
            // 相同则交换（gap 右移一格、当前列变成相等列），否则保持原样
            if query_base == b'-' {
                shift_gap_right(query_aligned, column, target_base);
            } else {
                shift_gap_right(target_aligned, column, query_base);
            }
            // 写入（可能已被交换过的）当前列；后续列会继续看到交换的结果
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

    /// 文档示例 (4) 的回归锚：mismatch 预支 t[2] 的 G 之后，列 2 的单侧 gap
    /// 又交换 t[3] 的 G 上来，输出列 2、3 均为 (G,G)（K=2 下可锚定 backbone）。
    #[test]
    fn test_doc_example_lookahead_then_swap() {
        assert_eq!(normalize("AGGT", "ACGG"), ("A-GGT".into(), "ACGG-".into()));
    }

    /// 文档示例补充：交换会跳过连续 gap 找到可配对的碱基，gap 落到远处。
    #[test]
    fn test_doc_example_swap_over_gap_run() {
        // 列1 的 q gap 与 q[3]='T' 交换 → 列1 变 (T,T)；列2、3 成 (-,T)
        assert_eq!(normalize("A--T", "ATTT"), ("AT--".into(), "ATTT".into()));
    }
}
