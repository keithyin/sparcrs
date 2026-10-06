# Sparc (pure Rust)

A sparsity-based consensus algorithm for long erroneous sequencing reads,
reimplemented in 100% safe Rust with zero runtime dependencies.

This is a pure Rust rewrite of the original Sparc C++ implementation
(Ye C, Ma Z. 2016, "Sparc: a sparsity-based consensus algorithm for long
erroneous sequencing reads", PeerJ 4:e2016, https://doi.org/10.7717/peerj.2016),
which previously lived at https://github.com/keithyin/Sparc as a C++ core with
Rust FFI bindings. The public API of the old binding crate is preserved.

- No C/C++ toolchain required — `cargo build` is all you need.
- No `unsafe`, no FFI, no external dependencies.
- Behavior-equivalent to the C++ implementation (validated by differential
  testing against the original binding); see "Behavior parity" below.
- 1.4–2.1x faster than the C++ binding on benchmark workloads.

## Usage

```rust
use sparc::{parse_m5, sparc_consensus, Query, SparcConfig};
use std::io::BufReader;

// backbone: reference-like sequence to polish
let backbone = "GATCGGGCTAA";

// each Query is one read alignment against the backbone
// (m5-style aligned strings with '-' gaps)
let queries = vec![
    Query::new(backbone.to_string(), "GATCGCGCTAA".to_string(), 0, backbone.len()),
    Query::new(backbone.to_string(), "GATCGCGCCAA".to_string(), 0, backbone.len()),
    Query::new(backbone.to_string(), "GCTCGGCCCAA".to_string(), 0, backbone.len()),
];

let mut config = SparcConfig::default();
config.kmer = 1;                // k-mer size, 1..=16 (suggested [1, 2])
config.coverage_threshold = 2;  // CLI "c", suggested [1, 5]
config.threshold = -0.1;        // CLI "t", adaptive threshold (<0 disables)

let consensus = sparc_consensus(backbone, &queries, &config).unwrap();
println!("consensus: {}", consensus.seq);
// consensus.start / consensus.end: best-path range on the backbone [start, end);
// end == 0 means "no trusted path, seq is the backbone as-is";
// start == None means the path head is not on the backbone.

// alignments can also be parsed from blasr m5 rows/files:
let m5 = std::fs::File::open("backbone-0.mapped.m5").unwrap();
let queries = parse_m5(BufReader::new(m5)).unwrap();
// negative-strand rows (tStrand '-') are handled automatically
```

All inputs are validated before the algorithm runs; invalid input (k-mer
out of range, backbone shorter than k, non-ACGT bases, inconsistent alignment
lengths/spans, out-of-range coordinates) returns a `SparcError` instead of
misbehaving.

## Parameters (CLI heritage)

| parameter | field | meaning |
|-----------|-------|---------|
| `k` | `kmer` | k-mer size (suggested [1, 2]) |
| `c` | `coverage_threshold` | coverage threshold (range [1, 5], suggest 2) |
| `t` | `threshold` | adaptive threshold ([0.0, 0.3]); `< 0` disables adaptivity |
| — | `cov_radius` | coverage sliding-window radius (original CLI: 200) |
| — | `scoring_method` | `ScoringMethod::Linear` (default) or `ScoringMethod::LogRatio` (reproduces a legacy C++ behavior) |
| — | `debug` | write debug artifacts into the process CWD |

## Performance

Measured with criterion (`cargo bench -p parity_tests`), release profile
(thin LTO, `codegen-units = 1`) vs. the C++ binding built by its `make -O3`:

| workload | pure Rust | C++ binding | speedup |
|----------|-----------|-------------|---------|
| 69 bp backbone, 7 reads (testdata2-like) | 35 µs | 50 µs | 1.4x |
| 10 kb backbone, 1000 × 300 bp reads (~30x) | 62 ms | 112 ms | 1.8x |
| 10 kb backbone, 5000 × 200 bp reads (~100x) | 280 ms | 575 ms | 2.1x |

The wins come from an arena-based graph (no per-node `malloc`), rolling
k-mer encoding (no per-k-mer bit-array surgery), and removing the FFI
boundary entirely.

## Behavior parity

The port is validated against the C++ binding by differential testing:
3860+ randomized scenarios (random backbones, substitutions/insertions/
deletions, reverse-strand alignments, orphan chains, k-mer collisions,
boundary radii/thresholds) must produce byte-identical
`Consensus { seq, start, end }` on both implementations
(`cargo test -p parity_tests`).

One deliberate deviation: the C++ core has a reachable out-of-bounds read
(`node_vec[MatchPosition]` in `GraphConstruction.cpp`, triggerable by
realistic alignments whose normalization creates double-gap columns; ASAN
confirms a heap-buffer-overflow, the plain build segfaults). The pure Rust
implementation treats that path as a defined no-op. The differential test
harness runs the C++ oracle in a subprocess and records such crashes as
known C++ bugs rather than parity failures.

Also not ported (dead code in the C++ core): MurmurHash64A/B,
`SparseConsensus*` structures, `Read`/`reads_table`/`Contigs`,
`SparcMergeNodes`, `SparcMultiply` (the `ScoringMethod::LogRatio` branch is
preserved bug-compatible), and the unused `kmer_t*` variants.

## Repository layout

```
sparc/          the publishable crate (name: sparc)
parity_tests/   dev-only crate: differential tests + benchmarks vs the C++ binding
sparc/testdata2/  golden test data
```

`parity_tests` needs the original C++ binding checked out at `../Sparc`
as its differential oracle; the `sparc` crate itself is fully standalone.

### Development notes

- `cargo fmt --all` also formats local path dependencies (documented cargo
  behavior) and would touch the oracle repository. Use `cargo fmt -p sparc`
  / `-p parity_tests` instead.
- The differential tests spawn the C++ oracle per scenario; they also serve
  as a regression suite for the known C++ UB documented above.

## Test data

https://sourceforge.net/projects/sparc-consensus/files/testdata/

## License

MIT
