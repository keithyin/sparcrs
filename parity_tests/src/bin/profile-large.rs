//! 性能剖析入口（开发专用）：构建 large 基准场景并反复跑 consensus，
//! 供 `perf record` 采样。不参与发布。

use std::hint::black_box;

use parity_tests::{Random, generate_alignment_with_span, random_base};

fn build_scenario(
    backbone_length: usize,
    read_count: usize,
    read_span: usize,
    seed: u64,
) -> (String, Vec<sparc::Query>, sparc::SparcConfig) {
    let mut rng = Random::new(seed);
    let mut backbone_bytes = Vec::with_capacity(backbone_length);
    for _ in 0..backbone_length {
        backbone_bytes.push(random_base(&mut rng));
    }
    let backbone = String::from_utf8(backbone_bytes).unwrap();

    let mut queries = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let alignment = generate_alignment_with_span(&mut rng, backbone.as_bytes(), read_span);
        queries.push(sparc::Query::new(
            alignment.target_aligned,
            alignment.query_aligned,
            alignment.target_start,
            alignment.target_end,
        ));
    }

    let config = sparc::SparcConfig {
        debug_output: None,
        kmer: 2,
        coverage_threshold: 2,
        scoring_method: sparc::ScoringMethod::Linear,
        cov_radius: 200,
        threshold: -0.1,
    };
    (backbone, queries, config)
}

fn main() {
    let (backbone, queries, config) = build_scenario(10_000, 5_000, 200, 0xC0FF_EE00 + 5_000);
    let iterations = std::env::args()
        .nth(1)
        .and_then(|v| v.parse::<u32>().ok())
        .unwrap_or(20);

    // 预热一次并打印结果规模，确保被优化器保留
    let first = sparc::sparc_consensus(&backbone, &queries, &config).expect("consensus 失败");
    eprintln!(
        "consensus len={} start={:?} end={}",
        first.seq.len(),
        first.start,
        first.end
    );

    let start = std::time::Instant::now();
    for _ in 0..iterations {
        let consensus =
            black_box(sparc::sparc_consensus(black_box(&backbone), &queries, &config)).unwrap();
        black_box(consensus.seq.len());
    }
    let elapsed = start.elapsed();
    eprintln!(
        "{iterations} iterations in {elapsed:.3?} ({:.1} ms/iter)",
        elapsed.as_millis() as f64 / f64::from(iterations)
    );
}
