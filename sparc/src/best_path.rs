//! 滑动窗口最大覆盖度与最优路径查找
//! （对应 C++ `sparc.cpp` 的 cov_cnt 滑窗 + `GraphSimplification.cpp` 的
//! `SparcFindBestPath`/`SparcBFSFindBestPath`）。

use std::collections::{BTreeMap, VecDeque};

use crate::graph::KmerGraph;
use crate::{ScoringMethod, SparcConfig};

/// 计算每个 backbone 位置的最大滑窗覆盖度（C++ `cov_vec`）。
///
/// 与 C++ 逐行对应的增量式计算：半径先 clamp（由调用方完成），
/// 窗口语义为 `[i - radius + 1, i + radius]`（移除先于加入），
/// 计数归零即从有序结构中移除，空窗口时记 0。
pub(crate) fn compute_window_max_coverage(
    nodes: &[crate::graph::GraphNode],
    radius: usize,
) -> Vec<u32> {
    let node_count = nodes.len();
    let mut window_max_coverage = vec![0u32; node_count];
    // C++ 用 map<int,int>；cov 值域内 BTreeMap 语义完全一致
    let mut coverage_counts: BTreeMap<u32, i64> = BTreeMap::new();

    for node in nodes.iter().take(radius) {
        *coverage_counts.entry(node.coverage).or_insert(0) += 1;
    }

    for index in 0..node_count {
        if index >= radius {
            let leaving_coverage = nodes[index - radius].coverage;
            let count = coverage_counts.entry(leaving_coverage).or_insert(0);
            *count -= 1;
            if *count <= 0 {
                coverage_counts.remove(&leaving_coverage);
            }
        }
        if index + radius < node_count {
            *coverage_counts
                .entry(nodes[index + radius].coverage)
                .or_insert(0) += 1;
        }
        window_max_coverage[index] = coverage_counts.keys().next_back().copied().unwrap_or(0);
    }

    window_max_coverage
}

/// 最优路径查找（C++ `SparcFindBestPath`）。
///
/// 依次以每个 backbone 节点（除最后一个）为源做松弛；得分跨源累积，
/// 等价于对整张图反复松弛直到收敛。分支计数等死代码不移植。
pub(crate) fn find_best_path(
    graph: &mut KmerGraph,
    config: &SparcConfig,
    window_max_coverage: &[u32],
) {
    if graph.backbone_node_count == 0 {
        return;
    }
    for (source_index, &source_window_max_coverage) in window_max_coverage
        [..graph.backbone_node_count - 1]
        .iter()
        .enumerate()
    {
        relax_from_source(graph, source_index, config, source_window_max_coverage);
    }
}

/// 以单个 backbone 节点为源的 BFS 松弛（C++ `SparcBFSFindBestPath`）。
///
/// 无 visited 集合：目标节点得分提升即（重复）入队，得分严格递增且有界故终止。
fn relax_from_source(
    graph: &mut KmerGraph,
    source_index: usize,
    config: &SparcConfig,
    source_window_max_coverage: u32,
) {
    let mut queue: VecDeque<usize> = VecDeque::new();
    queue.push_back(source_index);

    // 自适应阈值只依赖源位置的窗口最大覆盖度，C++ 在每条边上重复计算，
    // 这里提前到循环外（数值一致）
    let adaptive_threshold: i32 = if config.threshold < 0.0 {
        0 // 不使用
    } else {
        let mut value =
            ((source_window_max_coverage as i32) as f64 * config.threshold).round() as i32;
        if value < config.coverage_threshold {
            value = config.coverage_threshold;
        }
        value
    };

    while let Some(current_index) = queue.pop_front() {
        // 松弛过程不增删节点/边，先取边数再按下标访问，
        // 避免在迭代边的同时可变借用 nodes
        let edge_count = graph.nodes[current_index].right_edges.len();
        for edge_position in 0..edge_count {
            let (target_index, edge_coverage) = {
                let edge = &graph.nodes[current_index].right_edges[edge_position];
                (edge.target_node, edge.coverage)
            };

            // C++ 的 new_score 是 int（32 位），赋值时从 int64 截断；
            // 实际得分远小于 2^31，截断仅作语义对齐
            let new_score: i64;
            let updated: bool;

            match config.scoring_method {
                ScoringMethod::LogRatio => {
                    // C++ 方法 1 的既有行为：用大数比较决定是否更新，
                    // 但 new_score 恒为 0 —— 更新会把得分清零。按原样保留。
                    let current_score = graph.nodes[current_index].score;
                    let current_coverage = graph.nodes[current_index].coverage;
                    let target_score = graph.nodes[target_index].score;
                    let target_coverage = graph.nodes[target_index].coverage;
                    if target_score == 0 {
                        updated = true;
                    } else {
                        // C++ 比较 score/target_cov 与 (score+edge_cov)/current_cov
                        // 的交叉乘积（十进制大数字符串比较 = 数值比较）
                        updated = (target_score as i128) * (current_coverage as i128 + 1)
                            < ((current_score + edge_coverage as i64) as i128)
                                * (target_coverage as i128);
                    }
                    new_score = 0;
                }
                ScoringMethod::Linear => {
                    let current_score = graph.nodes[current_index].score;
                    let target_score = graph.nodes[target_index].score;
                    let candidate = if config.threshold < 0.0 {
                        // 每条边惩罚下限 -2
                        (current_score
                            + ((edge_coverage as i64) - (config.coverage_threshold as i64)).max(-2))
                            as i32
                    } else {
                        (current_score + (edge_coverage as i64) - (adaptive_threshold as i64))
                            as i32
                    };
                    new_score = candidate as i64;
                    updated = target_score < new_score;
                }
            }

            if updated {
                let target = &mut graph.nodes[target_index];
                target.score = new_score;
                target.best_predecessor = Some(current_index);
                if !target.in_backbone {
                    queue.push_back(target_index);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{GraphEdge, GraphNode};

    fn node_with_coverage(coverage: u32) -> GraphNode {
        GraphNode {
            kmer: 0,
            backbone_position: 0,
            score: 0,
            coverage,
            in_backbone: true,
            selected: false,
            best_predecessor: None,
            right_edges: Vec::new(),
            left_edges: Vec::new(),
        }
    }

    #[test]
    fn test_window_max_coverage_matches_incremental_semantics() {
        // 窗口语义 [i-radius+1, i+radius]（移除先于加入）
        let nodes: Vec<GraphNode> = [1u32, 5, 2]
            .iter()
            .map(|&c| node_with_coverage(c))
            .collect();
        assert_eq!(compute_window_max_coverage(&nodes, 1), vec![5, 5, 2]);

        // radius = 0：每步"移除自身再加入自身"，旧条目永远留在计数表里，
        // 退化为前缀运行最大值（与 C++ 增量算法的既有行为一致）
        assert_eq!(compute_window_max_coverage(&nodes, 0), vec![1, 5, 5]);

        // radius 覆盖全部 → 全局最大
        assert_eq!(compute_window_max_coverage(&nodes, 3), vec![5, 5, 5]);
    }

    fn build_test_graph() -> KmerGraph {
        // backbone "ACGT" (K=2)：3 个节点，0→1→2 主链，0→2 额外一条高覆盖边
        let mut graph = KmerGraph::build_backbone_chain("ACGT", 2);
        graph.nodes[0].right_edges.push(GraphEdge {
            target_node: 2,
            coverage: 10,
        });
        graph
    }

    #[test]
    fn test_relaxation_fixed_threshold() {
        let mut graph = build_test_graph();
        let config = SparcConfig {
            coverage_threshold: 2,
            threshold: -0.1,
            scoring_method: ScoringMethod::Linear,
            ..SparcConfig::default()
        };
        find_best_path(&mut graph, &config, &[5, 5, 5]);

        // 源 0：node1 得 0 + max(1-2, -2) = -1（不更新，0 < -1 为假），
        // node2 得 0 + max(10-2, -2) = 8 → 更新
        assert_eq!(graph.nodes[1].score, 0);
        assert_eq!(graph.nodes[2].score, 8);
        assert_eq!(graph.nodes[2].best_predecessor, Some(0));
        // 源 1：node2 = 0 + max(1-2, -2) = -1，不如已有 8，不更新
        assert_eq!(graph.nodes[2].best_predecessor, Some(0));
    }

    #[test]
    fn test_relaxation_adaptive_threshold() {
        let mut graph = build_test_graph();
        let config = SparcConfig {
            coverage_threshold: 2,
            threshold: 0.2,
            scoring_method: ScoringMethod::Linear,
            ..SparcConfig::default()
        };
        // 窗口最大覆盖度 5 → 阈值 = round(5 * 0.2) = 1 < CovTh → 取 2
        find_best_path(&mut graph, &config, &[5, 5, 5]);
        // node2 = 0 + 10 - 2 = 8
        assert_eq!(graph.nodes[2].score, 8);
    }

    #[test]
    fn test_log_ratio_writes_zero_scores() {
        // C++ 方法 1（LogRatio）的既有行为：更新时 new_score 恒为 0
        let mut graph = build_test_graph();
        let config = SparcConfig {
            scoring_method: ScoringMethod::LogRatio,
            threshold: -0.1,
            ..SparcConfig::default()
        };
        find_best_path(&mut graph, &config, &[5, 5, 5]);
        // target score == 0 时强制更新 → score 仍为 0，但 best_predecessor 被设置；
        // 源 1 的松弛会把 node2 的前驱覆盖为 1
        assert_eq!(graph.nodes[1].score, 0);
        assert_eq!(graph.nodes[1].best_predecessor, Some(0));
        assert_eq!(graph.nodes[2].score, 0);
        assert_eq!(graph.nodes[2].best_predecessor, Some(1));
    }
}
