//! 性能剖析入口（开发专用）：构建 large 基准场景并反复跑 consensus，
//! 供 `perf record` 采样。不参与发布。
//!
//! 场景定义与基准共用（`build_bench_scenario`），release 配置也与根
//! workspace 对齐（见 `parity_tests/Cargo.toml`）。用法：
//!
//! ```text
//! cargo run --manifest-path parity_tests/Cargo.toml --release \
//!     --bin profile-large -- [迭代次数]
//! ```

use std::hint::black_box;

use parity_tests::{BuiltScenario, LARGE_SCENARIO, build_bench_scenario};

fn main() {
    if cfg!(debug_assertions) {
        eprintln!("警告：当前为 debug 构建，剖析结果不代表 release 性能；请加 --release。");
    }

    // 用 args_os 避免非 UTF-8 参数触发 panic；非法输入直接报错而非静默取默认值。
    let iterations: u32 = match std::env::args_os().nth(1) {
        None => 20,
        Some(argument) => match argument.to_str().and_then(|value| value.parse().ok()) {
            Some(count) if count > 0 => count,
            _ => {
                eprintln!("用法: profile-large [迭代次数]（正整数，默认 20）");
                std::process::exit(2);
            }
        },
    };

    let BuiltScenario {
        backbone,
        config,
        queries_new: queries,
        ..
    } = build_bench_scenario(LARGE_SCENARIO);

    // 预热一次并打印结果规模，确保被优化器保留
    let first = sparc::sparc_consensus(&backbone, &queries, &config).expect("consensus 失败");
    eprintln!(
        "{} consensus len={} start={:?} end={}",
        LARGE_SCENARIO.name,
        first.seq.len(),
        first.start,
        first.end
    );

    let start = std::time::Instant::now();
    for _ in 0..iterations {
        let consensus = black_box(sparc::sparc_consensus(
            black_box(&backbone),
            black_box(queries.as_slice()),
            black_box(&config),
        ))
        .expect("consensus 失败");
        black_box(consensus.seq.len());
    }
    let elapsed = start.elapsed();
    let ms_per_iter = elapsed.as_secs_f64() * 1000.0 / f64::from(iterations);
    eprintln!("{iterations} iterations in {elapsed:.3?} ({ms_per_iter:.1} ms/iter)");
}
