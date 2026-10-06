# Changelog

All notable changes to the `sparc` crate are documented here.

## 0.7.0 (2026-10-06)

Breaking changes — API types and validation only; the computation itself is
unchanged and remains byte-identical to the C++ implementation across the
3860+ differential parity scenarios.

### Changed

- `SparcConfig` field types are now precise, making invalid states
  unrepresentable: `kmer: u8`, `coverage_threshold: u32`, `cov_radius: usize`.
  `SparcError::InvalidCovRadius` is gone (negative radius can no longer be
  expressed).
- `threshold = +inf` is rejected (`SparcError::InvalidThreshold`, message now
  "must not be NaN or +infinity"); negative values including `-inf` keep the
  "disable adaptive threshold" semantics.
- `debug: bool` is replaced by `debug_output: Option<DebugConfig>`. Debug
  artifacts (`align_profile.txt`, `subgraph.dot`, `subgraph_cns.dot`,
  `DEBUG.consensus.fasta`) are written into `DebugConfig::output_dir` instead
  of the process CWD, and file-creation failures now return the new
  `SparcError::DebugWrite` instead of panicking. The `subgraph_begin` /
  `subgraph_end` parameters moved into `DebugConfig` (as `usize`).
- `SparcError::M5Parse(String)` became `M5Parse { line, message }` with a
  1-based line number; parse errors now tell you which line failed.
- `SparcConfig` is no longer `Copy` (it holds `Option<DebugConfig>`); it
  remains `Clone`.

### Added

- `M5Reader<R: BufRead>` — an `Iterator<Item = Result<Query, SparcError>>`
  for streaming m5 files; `parse_m5` is now a thin wrapper over it.
- `Query::new` accepts `impl Into<String>`; read-only accessors for all
  `Query` fields (`query_aligned_sequence()`, `target_start()`, ...).
  `query_start`/`query_end` are documented as informational only (the
  algorithm does not use them).
- `Consensus::region() -> PathRegion` (`Backbone { start, end }` /
  `OrphanHead { end }` / `Fallback`) replaces manual decoding of the
  `start`/`end` sentinel encoding (which is kept for compatibility).
- `serde` feature (off by default) implementing `Serialize`/`Deserialize`
  for the public types.
- `#[must_use]` on `Query::reverse_strand`; `#![warn(missing_docs)]`;
  doctests for `Query`, `M5Reader`.
- CI: MSRV job (1.85), rustdoc with `-D warnings`, cargo-semver-checks,
  macOS build, `--features serde` test run, clippy over all targets.

## 0.6.0

- Idiomatic API: `ScoringMethod` enum, `usize` node indices, dropped reserved
  config fields; thiserror-based `SparcError`. (See git history.)
