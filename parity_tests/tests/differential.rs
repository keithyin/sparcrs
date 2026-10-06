//! 差分对拍测试：纯 Rust 实现与 C++ binding 的行为等价性。
//!
//! 每个 scenario 通过子进程调用 C++ oracle（见 `src/bin/legacy-oracle.rs`）。
//! C++ 在部分合法输入上会因既有 UB（`node_vec[MatchPosition]` 越界，见
//! GraphConstruction.cpp:905）崩溃——oracle 子进程崩溃会被检测并记录为
//! `legacy_crashed`，不算对拍失败（Rust 侧对该路径有安全定义的行为）。

use parity_tests::{Random, Scenario, generate_scenario, parse_oracle_line, serialize_scenario};
use std::io::Write;
use std::process::{Command, Stdio};

/// 一次 oracle 调用的结果。
enum OracleOutcome {
    Match(Result<(String, Option<u32>, u32), String>),
    Crashed,
}

fn run_legacy_oracle(scenario: &Scenario) -> OracleOutcome {
    let mut child = Command::new(env!("CARGO_BIN_EXE_legacy-oracle"))
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("启动 legacy-oracle 失败");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(serialize_scenario(scenario).as_bytes())
        .expect("写入 oracle stdin 失败");
    let output = child.wait_with_output().expect("等待 oracle 失败");
    if !output.status.success() {
        return OracleOutcome::Crashed;
    }
    let line = String::from_utf8_lossy(&output.stdout);
    match parse_oracle_line(&line) {
        Ok(result) => OracleOutcome::Match(result),
        Err(message) => panic!("oracle 输出异常: {message}"),
    }
}

/// 对拍一个 scenario。返回 Err(差异描述)。
fn compare_scenario(scenario: &Scenario) -> Result<bool, String> {
    // 先跑新实现（进程内，panic 会让测试失败）
    let new_result =
        sparc::sparc_consensus(&scenario.backbone, &scenario.queries_new, &scenario.config);

    match run_legacy_oracle(scenario) {
        OracleOutcome::Crashed => {
            assert!(
                new_result.is_ok(),
                "C++ oracle 崩溃的场景上，Rust 实现也必须给出结果，实际: {:?}",
                new_result
            );
            Ok(true) // legacy 崩溃（C++ 既有 UB），不算差异
        }
        OracleOutcome::Match(Ok((legacy_seq, legacy_start, legacy_end))) => {
            let new_consensus = new_result.expect("legacy 成功而新实现失败");
            if new_consensus.seq != legacy_seq {
                return Err(format!(
                    "seq 不一致\n  new    = {}\n  legacy = {}",
                    new_consensus.seq, legacy_seq
                ));
            }
            if new_consensus.start != legacy_start {
                return Err(format!(
                    "start 不一致: new={:?} legacy={:?}",
                    new_consensus.start, legacy_start
                ));
            }
            if new_consensus.end != legacy_end {
                return Err(format!(
                    "end 不一致: new={} legacy={}",
                    new_consensus.end, legacy_end
                ));
            }
            Ok(false)
        }
        OracleOutcome::Match(Err(legacy_message)) => match new_result {
            Err(new_error) => {
                if new_error.to_string() != legacy_message {
                    return Err(format!(
                        "错误消息不一致\n  new    = {new_error}\n  legacy = {legacy_message}"
                    ));
                }
                Ok(false)
            }
            Ok(consensus) => Err(format!(
                "legacy 报错而新实现成功: legacy={legacy_message}, new=Ok(start={:?}, end={})",
                consensus.start, consensus.end
            )),
        },
    }
}

fn run_differential(name: &str, seed_count: u64, scenario_builder: impl Fn(u64) -> Scenario) {
    let mut failures: Vec<String> = Vec::new();
    let mut legacy_crashed: Vec<u64> = Vec::new();
    for seed in 0..seed_count {
        let scenario = scenario_builder(seed);
        match compare_scenario(&scenario) {
            Ok(true) => legacy_crashed.push(seed),
            Ok(false) => {}
            Err(message) => {
                failures.push(format!("seed {seed}: {message}"));
                if failures.len() >= 5 {
                    break;
                }
            }
        }
    }
    eprintln!(
        "[{name}] 共 {seed_count} 个场景，其中 legacy（C++）崩溃 {} 个: {:?}",
        legacy_crashed.len(),
        legacy_crashed
    );
    assert!(
        failures.is_empty(),
        "{name}: {} 个场景不一致:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// 主对拍：随机 backbone/reads/参数。
#[test]
fn differential_parity_random_scenarios() {
    run_differential("random", 3000, |seed| {
        let mut rng = Random::new(seed.wrapping_mul(0x2545_F491_4F6C_DD1D) ^ 0xDEAD_BEEF);
        generate_scenario(&mut rng)
    });
}

/// 小 backbone（长度 == kmer，node_vec 只有 1 个节点）与边界配置。
#[test]
fn differential_parity_edge_cases() {
    run_differential("edge", 800, |seed| {
        let mut rng = Random::new(seed ^ 0x5EED_5EED);
        let backbone_length = 1 + rng.below(8);
        let mut backbone_bytes = Vec::with_capacity(backbone_length);
        for _ in 0..backbone_length {
            backbone_bytes.push(parity_tests::random_base(&mut rng));
        }
        let backbone = String::from_utf8(backbone_bytes).unwrap();

        let read_count = rng.below(12);
        let config = parity_tests::generate_config(&mut rng, backbone.len());
        let mut scenario = Scenario {
            backbone: backbone.clone(),
            legacy_config: parity_tests::to_legacy_config(&config),
            config,
            queries_new: Vec::new(),
            queries_legacy: Vec::new(),
            raw_alignments: Vec::new(),
        };

        for _ in 0..read_count {
            let alignment = parity_tests::generate_alignment(&mut rng, backbone.as_bytes());
            let new_query = sparc::Query::new(
                alignment.target_aligned.clone(),
                alignment.query_aligned.clone(),
                alignment.target_start,
                alignment.target_end,
            );
            let new_query = if alignment.reverse_strand {
                new_query.reverse_strand()
            } else {
                new_query
            };
            let legacy_query = sparc_legacy::Query::new(
                alignment.target_aligned.clone(),
                alignment.query_aligned.clone(),
                alignment.target_start,
                alignment.target_end,
            );
            let legacy_query = if alignment.reverse_strand {
                legacy_query.reverse_strand()
            } else {
                legacy_query
            };
            scenario.queries_new.push(new_query);
            scenario.queries_legacy.push(legacy_query);
            scenario.raw_alignments.push(alignment);
        }
        scenario
    });
}

/// 大规模场景（长 backbone、高覆盖、大半径）抽样。
#[test]
fn differential_parity_large_scenarios() {
    run_differential("large", 60, |seed| {
        let mut rng = Random::new(seed.wrapping_mul(0x9E37_79B9) | 0x8000_0000);
        let backbone_length = 1000 + rng.below(2000);
        let mut backbone_bytes = Vec::with_capacity(backbone_length);
        for _ in 0..backbone_length {
            backbone_bytes.push(parity_tests::random_base(&mut rng));
        }
        let backbone = String::from_utf8(backbone_bytes).unwrap();

        let mut scenario = Scenario {
            backbone: backbone.clone(),
            config: sparc::SparcConfig {
                debug: false,
                kmer: rng.pick(&[1, 2, 3]),
                coverage_threshold: 2,
                scoring_method: 2,
                subgraph_begin: 0,
                subgraph_end: backbone.len() as i32,
                cns_start: 0,
                cns_end: backbone.len() as i32,
                report_begin: 0,
                report_end: backbone.len() as i32,
                cov_radius: 200,
                threshold: -0.1,
            },
            legacy_config: sparc_legacy::SparcConfig::default(),
            queries_new: Vec::new(),
            queries_legacy: Vec::new(),
            raw_alignments: Vec::new(),
        };
        scenario.legacy_config = parity_tests::to_legacy_config(&scenario.config);

        for _ in 0..200 {
            let alignment = parity_tests::generate_alignment(&mut rng, backbone.as_bytes());
            let new_query = sparc::Query::new(
                alignment.target_aligned.clone(),
                alignment.query_aligned.clone(),
                alignment.target_start,
                alignment.target_end,
            );
            let new_query = if alignment.reverse_strand {
                new_query.reverse_strand()
            } else {
                new_query
            };
            let legacy_query = sparc_legacy::Query::new(
                alignment.target_aligned.clone(),
                alignment.query_aligned.clone(),
                alignment.target_start,
                alignment.target_end,
            );
            let legacy_query = if alignment.reverse_strand {
                legacy_query.reverse_strand()
            } else {
                legacy_query
            };
            scenario.queries_new.push(new_query);
            scenario.queries_legacy.push(legacy_query);
            scenario.raw_alignments.push(alignment);
        }
        scenario
    });
}
