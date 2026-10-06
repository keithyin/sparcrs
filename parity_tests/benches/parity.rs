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
use parity_tests::{Random, generate_alignment_with_span, random_base};

/// 生成均匀铺满 backbone 的高覆盖 reads。
fn build_scenario(
    backbone_length: usize,
    read_count: usize,
    read_span: usize,
    seed: u64,
) -> (
    String,
    Vec<sparc::Query>,
    Vec<sparc_legacy::Query>,
    sparc::SparcConfig,
) {
    let mut rng = Random::new(seed);
    let mut backbone_bytes = Vec::with_capacity(backbone_length);
    for _ in 0..backbone_length {
        backbone_bytes.push(random_base(&mut rng));
    }
    let backbone = String::from_utf8(backbone_bytes).unwrap();

    let mut queries_new = Vec::with_capacity(read_count);
    let mut queries_legacy = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let alignment = generate_alignment_with_span(&mut rng, backbone.as_bytes(), read_span);
        let new_query = sparc::Query::new(
            alignment.target_aligned.clone(),
            alignment.query_aligned.clone(),
            alignment.target_start,
            alignment.target_end,
        );
        queries_new.push(new_query);
        let legacy_query = sparc_legacy::Query::new(
            alignment.target_aligned,
            alignment.query_aligned,
            alignment.target_start,
            alignment.target_end,
        );
        queries_legacy.push(legacy_query);
    }

    let config = sparc::SparcConfig {
        debug: false,
        kmer: 2,
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
    };
    (backbone, queries_new, queries_legacy, config)
}

fn bench_against_legacy(criterion: &mut Criterion) {
    let scenarios = [
        ("small_69bp_7reads", 69, 7, 60),
        ("medium_10kb_1000reads_30x", 10_000, 1_000, 300),
        ("large_10kb_5000reads_100x", 10_000, 5_000, 200),
    ];

    let mut group = criterion.benchmark_group("consensus");
    for (name, backbone_length, read_count, read_span) in scenarios {
        let (backbone, queries_new, queries_legacy, config) = build_scenario(
            backbone_length,
            read_count,
            read_span,
            0xC0FF_EE00 + read_count as u64,
        );
        let legacy_config = parity_tests::to_legacy_config(&config);

        // 基准中不做断言，但先各跑一次确保两边都正常出结果
        let new_result =
            sparc::sparc_consensus(&backbone, &queries_new, &config).expect("新实现失败");
        let legacy_result =
            sparc_legacy::sparc_consensus(&backbone, &queries_legacy, &legacy_config)
                .expect("C++ 实现失败");
        assert_eq!(new_result.seq, legacy_result.seq, "基准场景两边结果不一致");
        assert_eq!(new_result.start, legacy_result.start);
        assert_eq!(new_result.end, legacy_result.end);

        let sample_size = if read_count > 1_000 { 20 } else { 50 };
        group.sample_size(sample_size);

        group.bench_with_input(
            BenchmarkId::new("rust", name),
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
            BenchmarkId::new("cpp", name),
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
