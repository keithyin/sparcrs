//! 新旧实现的性能对比基准（三档负载）。
//!
//! 运行：`cargo bench -p parity_tests`
//!
//! | 场景  | backbone | reads          | 覆盖度 |
//! |-------|----------|----------------|--------|
//! | small | 69bp     | 7 × ~60bp      | ~6x    |
//! | medium| 10kb     | 1000 × 300bp   | ~30x   |
//! | large | 10kb     | 5000 × 200bp   | ~100x  |

use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use parity_tests::{BENCH_SCENARIOS, BuiltScenario, build_bench_scenario, to_legacy_config};

fn bench_against_legacy(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("consensus");
    for scenario in BENCH_SCENARIOS {
        let BuiltScenario {
            backbone,
            config,
            queries_new,
            queries_legacy,
        } = build_bench_scenario(scenario);
        let legacy_config = to_legacy_config(&config);

        // 基准中不做断言，但先各跑一次确保两边都正常出结果
        let new_result =
            sparc::sparc_consensus(&backbone, &queries_new, &config).expect("新实现失败");
        let legacy_result =
            sparc_legacy::sparc_consensus(&backbone, &queries_legacy, &legacy_config)
                .expect("C++ 实现失败");
        assert_eq!(new_result.seq, legacy_result.seq, "基准场景两边结果不一致");
        // legacy 的 start/end 为 u32（C++ uint32_t 语义），在对拍边界转换
        assert_eq!(new_result.start, legacy_result.start.map(|v| v as usize));
        assert_eq!(new_result.end, legacy_result.end as usize);

        let sample_size = if scenario.read_count > 1_000 { 20 } else { 50 };
        group.sample_size(sample_size);

        group.bench_with_input(
            BenchmarkId::new("rust", scenario.name),
            &(&backbone, &queries_new, &config),
            |benchmark, (backbone, queries, config)| {
                benchmark.iter(|| {
                    black_box(sparc::sparc_consensus(
                        black_box(backbone),
                        black_box(queries),
                        black_box(config),
                    ))
                    .expect("新实现失败")
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("cpp", scenario.name),
            &(&backbone, &queries_legacy, &legacy_config),
            |benchmark, (backbone, queries, config)| {
                benchmark.iter(|| {
                    black_box(sparc_legacy::sparc_consensus(
                        black_box(backbone),
                        black_box(queries),
                        black_box(config),
                    ))
                    .expect("C++ 实现失败")
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_against_legacy);
criterion_main!(benches);
