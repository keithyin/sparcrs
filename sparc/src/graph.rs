//! k-mer 图（对应 C++ `GraphConstruction.cpp` + `ConsensusNode`/`ConsensusEdgeNode`）。
//!
//! 内存模型：arena（`Vec<GraphNode>` + u32 索引）替代 C++ 的逐节点 malloc +
//! 裸指针链表。前 `backbone_node_count` 个节点即 C++ 的 `node_vec`
//! （按 backbone k-mer 起始位置索引），其后是分支/孤儿节点。
//! 孤儿节点在 C++ 中需要 `orphan_nodes` 单独登记才能安全释放；
//! arena 下生命周期由 Vec 统一管理，无需登记。
//!
//! 边匹配语义（必须与 C++ 完全一致，两种区域不同）：
//! - match 区：按**节点身份**匹配（C++ `node_ptr == current_node`）；
//! - 分支区：right 边按「目标节点 kmer 值 + 非 backbone」匹配，
//!   left 边按「目标节点 kmer 值」匹配（C++ GraphConstruction.cpp:1051 的
//!   既有行为，即使目标节点不是真正的 previous 也计入其边覆盖度）。

use crate::Query;
use crate::encoding::{base_to_two_bit_code, reverse_complement_in_place};
use crate::normalize::normalize_alignment;

/// 图节点（对应 C++ `ConsensusNode`，去掉位域与无用字段）。
#[derive(Debug, Clone)]
pub(crate) struct GraphNode {
    /// 2-bit 编码的 k-mer（首碱基在高位），K ≤ 16 时不超过 32 bit。
    pub(crate) kmer: u32,
    /// backbone 位置（C++ `coord`）：backbone 节点 = 自身下标，分支节点 = 0。
    pub(crate) backbone_position: u32,
    /// 最优路径动态规划的得分（C++ `score`）。
    pub(crate) score: i64,
    /// 覆盖度（C++ `cov`）：backbone 节点初值 1，每次命中 +1。
    pub(crate) coverage: u32,
    /// 是否为 backbone 主链节点。
    pub(crate) in_backbone: bool,
    /// debug 输出用：是否在选出的共识路径上。
    pub(crate) selected: bool,
    /// 最优路径回溯指针（C++ `last_node`）；None = 无前驱。
    pub(crate) best_predecessor: Option<u32>,
    /// 右向边（到后继节点）。
    pub(crate) right_edges: Vec<GraphEdge>,
    /// 左向边（到前驱节点）。
    pub(crate) left_edges: Vec<GraphEdge>,
}

/// 图边（对应 C++ `ConsensusEdgeNode`）。
///
/// C++ 中 right/left 是两条独立的边链表（各自拥有独立的边对象与覆盖度，
/// 通常成对创建但可能分叉），这里用两个独立的 Vec 保留该语义。
#[derive(Debug, Clone)]
pub(crate) struct GraphEdge {
    /// 边指向的节点下标（C++ `node_ptr`）。
    pub(crate) target_node: u32,
    /// 边覆盖度（C++ `edge_cov`）。
    pub(crate) coverage: i32,
}

/// 整张 k-mer 图。
pub(crate) struct KmerGraph {
    /// 所有节点：前 `backbone_node_count` 个为 backbone 主链节点。
    pub(crate) nodes: Vec<GraphNode>,
    /// backbone 主链节点数（C++ `node_vec.size()` = backbone 长度 - k + 1）。
    pub(crate) backbone_node_count: usize,
}

impl KmerGraph {
    /// 构建 backbone 主链（对应 C++ `SparcConsensusKmerGraphConstruction`）。
    ///
    /// 每个 backbone 节点 coverage 初值 1、in_backbone = true、
    /// backbone_position = 自身下标；相邻节点以 coverage = 1 的双向边相连。
    pub(crate) fn build_backbone_chain(backbone: &str, kmer_size: usize) -> KmerGraph {
        let backbone_node_count = backbone.len() - kmer_size + 1;
        let mut nodes: Vec<GraphNode> = Vec::with_capacity(backbone_node_count);
        // 滚动 k-mer：等价于 C++ 逐位置 get_sub_arr 抽取 + str2bitsarr 编码
        let kmer_mask: u64 = (1u64 << (2 * kmer_size)) - 1;
        let mut rolling_kmer: u64 = 0;

        for (position, &base) in backbone.as_bytes().iter().enumerate() {
            rolling_kmer = ((rolling_kmer << 2) | base_to_two_bit_code(base)) & kmer_mask;
            if position + 1 < kmer_size {
                continue;
            }
            let node_index = position + 1 - kmer_size;
            let mut node = GraphNode {
                kmer: rolling_kmer as u32,
                backbone_position: node_index as u32,
                score: 0,
                coverage: 1,
                in_backbone: true,
                selected: false,
                best_predecessor: None,
                right_edges: Vec::new(),
                left_edges: Vec::new(),
            };
            if node_index >= 1 {
                nodes[node_index - 1].right_edges.push(GraphEdge {
                    target_node: node_index as u32,
                    coverage: 1,
                });
                node.left_edges.push(GraphEdge {
                    target_node: (node_index - 1) as u32,
                    coverage: 1,
                });
            }
            nodes.push(node);
        }

        KmerGraph {
            nodes,
            backbone_node_count,
        }
    }

    fn push_node(&mut self, kmer: u32, in_backbone: bool) -> u32 {
        assert!(
            self.nodes.len() < u32::MAX as usize,
            "k-mer 图节点数超出 u32 索引范围"
        );
        let node_index = self.nodes.len() as u32;
        self.nodes.push(GraphNode {
            kmer,
            backbone_position: 0,
            score: 0,
            coverage: 1,
            in_backbone,
            selected: false,
            best_predecessor: None,
            right_edges: Vec::new(),
            left_edges: Vec::new(),
        });
        node_index
    }

    /// 将一条 read 的比对路径加入图中（对应 C++ `SparcAddPathToBackbone`）。
    ///
    /// 库路径下 Patch/Fill 恒为 false（对应 C++ NewQuery 的默认值），
    /// 因此只做反向互补（负链）与归一化。
    pub(crate) fn add_path_to_backbone(&mut self, query: &Query, kmer_size: usize) {
        let mut query_aligned = query.query_aligned_sequence.bytes().collect::<Vec<u8>>();
        let mut target_aligned = query.target_aligned_sequence.bytes().collect::<Vec<u8>>();

        // 负链：两条比对串按 read 方向给出，先反向互补（C++ tStrand == '-' 分支）
        if query.is_reverse_strand {
            reverse_complement_in_place(&mut query_aligned);
            reverse_complement_in_place(&mut target_aligned);
        }

        normalize_alignment(&mut query_aligned, &mut target_aligned);

        let query_length = query_aligned.len();
        let mut target_position: usize = query.target_start;
        let mut previous_node: Option<u32> = None;
        let mut current_node: Option<u32>;
        // C++ 的 MatchPosition：Some(backbone 位置) = match 区，None = 分支区
        let mut match_position: Option<usize>;

        let mut column = 0usize;
        while column + kmer_size <= query_length {
            // query 侧 gap：不产生节点，仅推进 target 位置。
            // （C++ 在此把 MatchPosition 置 -1，但下一轮循环在任何读取之前
            // 都会重新赋值，故该赋值是死代码，Rust 侧省略。）
            if query_aligned[column] == b'-' {
                if target_aligned[column] != b'-' {
                    target_position += 1;
                }
                column += 1;
                continue;
            }

            // 判断当前列开始的 K 个列是否全部相等（双 gap 列算相等）
            let mut window_is_match = true;
            for offset in 0..kmer_size {
                if query_aligned[column + offset] != target_aligned[column + offset] {
                    window_is_match = false;
                    break;
                }
            }
            match_position = if window_is_match {
                Some(target_position)
            } else {
                None
            };

            // 从当前列起收集 query 侧 K 个非 gap 碱基并编码为 k-mer
            let mut kmer: u64 = 0;
            let mut collected = 0usize;
            for &base in query_aligned[column..].iter() {
                if collected == kmer_size {
                    break;
                }
                if base != b'-' {
                    kmer = (kmer << 2) | base_to_two_bit_code(base);
                    collected += 1;
                }
            }
            if collected < kmer_size {
                // 串尾不足一个完整 k-mer，结束
                break;
            }
            let kmer = kmer as u32;

            match match_position {
                Some(position) => {
                    // match 区：锚定到 backbone 节点
                    // C++ 直接下标访问 node_vec[MatchPosition]；越界属于
                    // 病态输入（历史 UB 路径），这里取 None 跳过以保证安全
                    let anchored = self.nodes.get(position).map(|_| position as u32);
                    if let Some(anchored_index) = anchored {
                        self.nodes[anchored_index as usize].coverage += 1;
                        if let Some(previous_index) = previous_node {
                            upsert_edge_by_identity(
                                &mut self.nodes,
                                previous_index,
                                anchored_index,
                                Side::Right,
                            );
                            upsert_edge_by_identity(
                                &mut self.nodes,
                                anchored_index,
                                previous_index,
                                Side::Left,
                            );
                        }
                    }
                    current_node = anchored;
                }
                None => {
                    // 分支区（insertion / mismatch）
                    match previous_node {
                        None => {
                            // 链首孤儿节点：无入向边（C++ 登记进 orphan_nodes）
                            let orphan_index = self.push_node(kmer, false);
                            current_node = Some(orphan_index);
                        }
                        Some(previous_index) => {
                            // 在 previous 的右边里按「目标 kmer 值 + 非 backbone」找
                            let mut hit: Option<(usize, u32)> = None;
                            for (edge_position, edge) in self.nodes[previous_index as usize]
                                .right_edges
                                .iter()
                                .enumerate()
                            {
                                let target = &self.nodes[edge.target_node as usize];
                                if target.kmer == kmer && !target.in_backbone {
                                    hit = Some((edge_position, edge.target_node));
                                    break;
                                }
                            }

                            match hit {
                                Some((edge_position, hit_index)) => {
                                    self.nodes[previous_index as usize].right_edges
                                        [edge_position]
                                        .coverage += 1;
                                    self.nodes[hit_index as usize].coverage += 1;
                                    current_node = Some(hit_index);
                                }
                                None => {
                                    let branch_index = self.push_node(kmer, false);
                                    self.nodes[previous_index as usize].right_edges.push(
                                        GraphEdge {
                                            target_node: branch_index,
                                            coverage: 1,
                                        },
                                    );
                                    current_node = Some(branch_index);
                                }
                            }

                            // 回连 previous：按「目标节点 kmer 值 == previous 的 kmer 值」
                            // 匹配（C++ 原始行为，可能命中同 kmer 的其他节点）
                            let current_index = current_node.expect("current_node 已设置");
                            let previous_kmer = self.nodes[previous_index as usize].kmer;
                            let left_hit = self.nodes[current_index as usize]
                                .left_edges
                                .iter()
                                .position(|edge| {
                                    self.nodes[edge.target_node as usize].kmer == previous_kmer
                                });
                            match left_hit {
                                Some(edge_position) => {
                                    self.nodes[current_index as usize].left_edges[edge_position]
                                        .coverage += 1;
                                }
                                None => {
                                    self.nodes[current_index as usize]
                                        .left_edges
                                        .push(GraphEdge {
                                            target_node: previous_index,
                                            coverage: 1,
                                        });
                                }
                            }
                        }
                    }
                }
            }

            previous_node = current_node;
            if target_aligned[column] != b'-' {
                target_position += 1;
            }
            column += 1;
        }
    }
}

/// 边的挂载侧。
#[derive(Clone, Copy)]
enum Side {
    Right,
    Left,
}

/// 按「目标节点身份」查找边，命中则覆盖度 +1，否则新建边（coverage = 1）。
/// 对应 C++ match 区的边链表扫描（`node_ptr == current_node`）。
fn upsert_edge_by_identity(nodes: &mut [GraphNode], from_node: u32, to_node: u32, side: Side) {
    let edges = match side {
        Side::Right => &mut nodes[from_node as usize].right_edges,
        Side::Left => &mut nodes[from_node as usize].left_edges,
    };
    for edge in edges.iter_mut() {
        if edge.target_node == to_node {
            edge.coverage += 1;
            return;
        }
    }
    edges.push(GraphEdge {
        target_node: to_node,
        coverage: 1,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn query_from(query: &str, target: &str, target_start: usize, target_end: usize) -> Query {
        Query::new(
            target.to_string(),
            query.to_string(),
            target_start,
            target_end,
        )
    }

    #[test]
    fn test_build_backbone_chain() {
        let graph = KmerGraph::build_backbone_chain("ACGT", 2);
        assert_eq!(graph.backbone_node_count, 3);
        // kmer 值：AC = (0<<2)|1 = 1，CG = (1<<2)|2 = 6，GT = (2<<2)|3 = 11
        assert_eq!(graph.nodes[0].kmer, 1);
        assert_eq!(graph.nodes[1].kmer, 6);
        assert_eq!(graph.nodes[2].kmer, 11);
        for node in &graph.nodes {
            assert!(node.in_backbone);
            assert_eq!(node.coverage, 1);
            assert_eq!(node.score, 0);
        }
        // 链式边
        assert_eq!(graph.nodes[0].right_edges.len(), 1);
        assert_eq!(graph.nodes[0].right_edges[0].target_node, 1);
        assert_eq!(graph.nodes[0].right_edges[0].coverage, 1);
        assert_eq!(graph.nodes[1].left_edges.len(), 1);
        assert_eq!(graph.nodes[1].left_edges[0].target_node, 0);
        assert_eq!(graph.nodes[2].left_edges[0].target_node, 1);
        // backbone_position = 自身下标
        assert_eq!(graph.nodes[0].backbone_position, 0);
        assert_eq!(graph.nodes[1].backbone_position, 1);
    }

    #[test]
    fn test_exact_match_query_doubles_coverage() {
        let mut graph = KmerGraph::build_backbone_chain("ACGT", 2);
        let query = query_from("ACGT", "ACGT", 0, 4);
        graph.add_path_to_backbone(&query, 2);
        // 每个节点被命中一次：coverage 1 → 2
        assert_eq!(graph.nodes[0].coverage, 2);
        assert_eq!(graph.nodes[1].coverage, 2);
        assert_eq!(graph.nodes[2].coverage, 2);
        // 边去重：不新增边，覆盖度翻倍
        assert_eq!(graph.nodes[0].right_edges.len(), 1);
        assert_eq!(graph.nodes[0].right_edges[0].coverage, 2);
        assert_eq!(graph.nodes[1].left_edges[0].coverage, 2);
        assert_eq!(graph.nodes.len(), 3);
    }

    #[test]
    fn test_insertion_creates_branch_node() {
        let mut graph = KmerGraph::build_backbone_chain("ACGT", 2);
        // 在 C 后插入 G：q = ACGGT, t = AC-GT（span 0..4）
        let query = query_from("ACGGT", "AC-GT", 0, 4);
        graph.add_path_to_backbone(&query, 2);
        // 归一化后 mismatch/gap 列调整，但应产生一个非 backbone 分支节点
        let branch_nodes: Vec<&GraphNode> =
            graph.nodes[graph.backbone_node_count..].iter().collect();
        assert!(!branch_nodes.is_empty());
        for node in &branch_nodes {
            assert!(!node.in_backbone);
            assert_eq!(node.backbone_position, 0);
        }
    }

    #[test]
    fn test_leading_mismatch_creates_orphan() {
        let mut graph = KmerGraph::build_backbone_chain("ACGT", 2);
        // 首碱基 mismatch：链首进入分支区且 previous == None → 孤儿节点
        let query = query_from("TCGT", "ACGT", 0, 4);
        graph.add_path_to_backbone(&query, 2);
        assert!(graph.nodes.len() > graph.backbone_node_count);
        let orphan = &graph.nodes[graph.backbone_node_count];
        assert!(!orphan.in_backbone);
        // 孤儿链首没有入向右边的来源，left 边为空
        assert!(orphan.left_edges.is_empty());
    }

    #[test]
    fn test_negative_strand_matches_forward() {
        let mut forward = KmerGraph::build_backbone_chain("ACGT", 2);
        let forward_query = query_from("ACGT", "ACGT", 0, 4);
        forward.add_path_to_backbone(&forward_query, 2);

        let mut reversed = KmerGraph::build_backbone_chain("ACGT", 2);
        // read 方向的负链比对串 = 正链串的反向互补
        let mut revcomp_q = b"ACGT".to_vec();
        let mut revcomp_t = b"ACGT".to_vec();
        reverse_complement_in_place(&mut revcomp_q);
        reverse_complement_in_place(&mut revcomp_t);
        let reversed_query = Query {
            query_aligned_sequence: String::from_utf8(revcomp_q).unwrap(),
            target_aligned_sequence: String::from_utf8(revcomp_t).unwrap(),
            is_reverse_strand: true,
            query_start: 0,
            query_end: 4,
            target_start: 0,
            target_end: 4,
        };
        reversed.add_path_to_backbone(&reversed_query, 2);

        let forward_covs: Vec<u32> = forward.nodes.iter().map(|n| n.coverage).collect();
        let reversed_covs: Vec<u32> = reversed.nodes.iter().map(|n| n.coverage).collect();
        assert_eq!(forward_covs, reversed_covs);
    }
}
