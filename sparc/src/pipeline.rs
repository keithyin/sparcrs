//! 端到端流水线编排：建图 → 加读段路径 → 滑窗覆盖度 → 最优路径 →
//! 选终点 → 回溯解码 → 输出（对应 C++ `sparc.cpp:SparcConsensus`）。
//!
//! 本模块只做编排，不做具体计算；各阶段分别位于
//! [`crate::graph`]、[`crate::best_path`] 等模块。

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::best_path::{compute_window_max_coverage, find_best_path};
use crate::encoding::decode_kmer_to_string;
use crate::graph::KmerGraph;
use crate::{Consensus, Query, SparcConfig};

/// 前置条件：输入已由 `validate_inputs` 校验通过。
pub(crate) fn compute_consensus(
    backbone: &str,
    queries: &[Query],
    config: &SparcConfig,
) -> Consensus {
    let kmer_size = config.kmer as usize;

    // 1. backbone 主链
    let mut graph = KmerGraph::build_backbone_chain(backbone, kmer_size);

    // debug 模式：align_profile.txt 在 C++ 中于建图前创建（可能最终为空文件）
    let mut align_profile = if config.debug {
        Some(BufWriter::new(
            File::create("align_profile.txt").expect("创建 align_profile.txt 失败"),
        ))
    } else {
        None
    };

    // 2. 逐条读段加入分支路径
    for query in queries {
        graph.add_path_to_backbone(query, kmer_size);
    }

    // 3. 滑动窗口最大覆盖度（半径两段 clamp，与 C++ 一致）
    let mut radius = config.cov_radius as usize;
    if radius > backbone.len() {
        radius = backbone.len() - 1;
    }
    if radius > graph.backbone_node_count {
        radius = graph.backbone_node_count;
    }
    let window_max_coverage = compute_window_max_coverage(&graph.nodes, radius);

    // 4. 最优路径
    find_best_path(&mut graph, config, &window_max_coverage);

    // debug：subgraph.dot 在找路之后、终点选择之前输出
    if config.debug && config.subgraph_end > config.subgraph_begin {
        write_subgraph_dot(
            &graph,
            config.subgraph_begin,
            config.subgraph_end,
            "subgraph.dot",
        );
    }

    // 5. 终点选择：backbone 节点中第一个得分最大者（严格大于）
    let backbone_node_count = graph.backbone_node_count;
    let mut max_score: i64 = 0;
    let mut end_position_inclusive: usize = 0;
    for index in 0..backbone_node_count {
        if graph.nodes[index].score > max_score {
            max_score = graph.nodes[index].score;
            end_position_inclusive = index;
        }
    }

    // fallback：无可信路径，原样返回 backbone
    if max_score == 0 {
        return Consensus {
            seq: backbone.to_string(),
            start: None,
            end: 0,
        };
    }

    // 6. 回溯解码（C++ 逐行对应）：
    //    末节点解码整串后先反转，前驱节点各贡献解码串首字符，最后整体反转
    let mut consensus_bytes: Vec<u8> = Vec::new();
    let mut start_node: Option<u32> = None;
    {
        let end_node = &graph.nodes[end_position_inclusive];
        let decoded = decode_kmer_to_string(end_node.kmer, kmer_size);
        consensus_bytes.extend(decoded.iter().rev());
        let mut cursor = end_node.best_predecessor;
        while let Some(node_index) = cursor {
            let node = &graph.nodes[node_index as usize];
            let decoded = decode_kmer_to_string(node.kmer, kmer_size);
            consensus_bytes.push(decoded[0]);
            if node.best_predecessor.is_none() {
                start_node = Some(node_index);
            }
            cursor = node.best_predecessor;
        }
    }
    consensus_bytes.reverse();

    // 路径头是 backbone 节点时其下标即 start，否则（孤儿链首）为 None
    let start_position =
        start_node.filter(|&node_index| (node_index as usize) < backbone_node_count);

    // 7. debug 输出
    if config.debug {
        write_debug_consensus(&consensus_bytes);
        mark_selected_path(&mut graph, end_position_inclusive);
        if config.subgraph_end > config.subgraph_begin {
            write_subgraph_dot(
                &graph,
                config.subgraph_begin,
                config.subgraph_end,
                "subgraph_cns.dot",
            );
        }
        if let Some(profile) = align_profile.as_mut() {
            for index in 0..backbone_node_count {
                let node = &graph.nodes[index];
                let _ = write!(profile, "{index} {}", node.score);
                if node.backbone_position > 0 {
                    let _ = write!(profile, " bb_coord: {}", node.backbone_position);
                }
                let _ = writeln!(profile, " cns_coord: 0");
            }
        }
    }

    Consensus {
        seq: String::from_utf8(consensus_bytes).expect("k-mer 解码结果必为 ASCII"),
        start: start_position,
        end: end_position_inclusive as u32 + 1,
    }
}

/// debug 模式：沿 best_predecessor 标记共识路径上的节点。
fn mark_selected_path(graph: &mut KmerGraph, end_position_inclusive: usize) {
    let mut cursor = Some(end_position_inclusive as u32);
    while let Some(node_index) = cursor {
        let node = &mut graph.nodes[node_index as usize];
        node.selected = true;
        cursor = node.best_predecessor;
    }
}

/// debug 模式：输出 DEBUG.consensus.fasta（格式与 C++ 一致）。
fn write_debug_consensus(consensus_bytes: &[u8]) {
    let consensus = String::from_utf8_lossy(consensus_bytes);
    if let Ok(output) = File::create("DEBUG.consensus.fasta") {
        let mut output = BufWriter::new(output);
        let _ = writeln!(output, ">Debug");
        let _ = writeln!(output, "{consensus}");
    }
}

/// debug 模式：输出子图 dot 文件（C++ 用指针地址作节点 id，这里用节点
/// 下标 `n{index}`，拓扑内容一致；越界区间按 C++ 的 UB 行为裁剪到有效范围）。
fn write_subgraph_dot(graph: &KmerGraph, subgraph_begin: i32, subgraph_end: i32, filename: &str) {
    let node_count = graph.nodes.len() as i32;
    let begin = subgraph_begin.clamp(0, node_count);
    let end = subgraph_end.clamp(0, node_count);
    if begin >= end {
        return;
    }
    let output = match File::create(filename) {
        Ok(output) => output,
        Err(_) => return,
    };
    let mut output = BufWriter::new(output);
    let _ = writeln!(output, "digraph G {{");

    for index in begin..end {
        let node = &graph.nodes[index as usize];
        for edge in &node.right_edges {
            let target = &graph.nodes[edge.target_node as usize];
            let _ = writeln!(
                output,
                "\"n{index}\" -> \"n{}\" [label=\"{}\"];",
                edge.target_node, edge.coverage
            );
            let base = crate::encoding::two_bit_code_to_base((target.kmer & 3) as u64) as char;
            let backbone_flag = if target.in_backbone { ":B" } else { "" };
            let selected_flag = if target.selected { "color=red," } else { "" };
            let _ = writeln!(
                output,
                "\"n{}\" [{selected_flag}label=\"{base}[{}][{}]{}\"];",
                edge.target_node, target.backbone_position, target.coverage, backbone_flag,
            );
        }
    }
    let _ = writeln!(output, "}}");
}
