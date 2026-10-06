//! legacy oracle：独立进程运行 C++ binding。
//!
//! 差分对拍时，C++ 实现在部分合法输入上会因既有 UB（GraphConstruction.cpp:905
//! 的 `node_vec[MatchPosition]` 越界）直接崩溃。把 oracle 放进子进程后，
//! 崩溃可以被父进程检测并记录，而不会杀死测试进程。
//!
//! 协议：stdin 读入一行场景（见 `parity_tests::serialize_scenario`），
//! stdout 输出一行 `OK \t start \t end \t seq` 或 `ERR \t message`。

use std::io::{BufRead, Write};

fn main() {
    let mut line = String::new();
    std::io::stdin()
        .lock()
        .read_line(&mut line)
        .expect("oracle 读取 stdin 失败");

    let fields: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
    let mut cursor = 0usize;
    let mut next = || {
        let value = fields[cursor];
        cursor += 1;
        value
    };

    let backbone = next().to_string();
    let config = sparc_legacy::SparcConfig {
        debug: false,
        kmer: next().parse().expect("kmer"),
        coverage_threshold: next().parse().expect("coverage_threshold"),
        scoring_method: next().parse().expect("scoring_method"),
        subgraph_begin: next().parse().expect("subgraph_begin"),
        subgraph_end: next().parse().expect("subgraph_end"),
        // 保留字段仅对原 CLI 读入的 m5 行生效，对 FFI/API 传入的 query 无影响
        cns_start: 0,
        cns_end: 0,
        report_begin: 0,
        report_end: 0,
        cov_radius: next().parse().expect("cov_radius"),
        threshold: next().parse().expect("threshold"),
    };
    let query_count: usize = next().parse().expect("query_count");

    let mut queries = Vec::with_capacity(query_count);
    for _ in 0..query_count {
        let target_start: usize = next().parse().expect("target_start");
        let target_end: usize = next().parse().expect("target_end");
        let reverse_strand = next() == "1";
        let query_aligned = next().to_string();
        let target_aligned = next().to_string();
        let query =
            sparc_legacy::Query::new(target_aligned, query_aligned, target_start, target_end);
        queries.push(if reverse_strand {
            query.reverse_strand()
        } else {
            query
        });
    }

    let result = sparc_legacy::sparc_consensus(&backbone, &queries, &config);
    let stdout = std::io::stdout();
    let mut stdout = stdout.lock();
    match result {
        Ok(consensus) => {
            let start = consensus
                .start
                .map(|value| value.to_string())
                .unwrap_or_else(|| "-1".into());
            writeln!(stdout, "OK\t{start}\t{}\t{}", consensus.end, consensus.seq).unwrap();
        }
        Err(error) => {
            writeln!(stdout, "ERR\t{error}").unwrap();
        }
    }
}
