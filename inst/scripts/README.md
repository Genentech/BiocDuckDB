# BiocDuckDB benchmark scripts

These scripts produce the numbers rendered in the *Benchmarking BiocDuckDB*
vignette, following the
[HDF5Array performance vignette](https://bioconductor.org/packages/release/bioc/vignettes/HDF5Array/inst/doc/HDF5Array_performance.html)
precompute pattern: the full-scale benchmark runs **offline** and writes a results
file that the vignette renders, so the vignette itself builds quickly.

## Files

- **`run_scran_scuttle_benchmarks.R`** — times the SQL-optimized `scuttle`/`scran`
  methods (QC, normalization, pseudo-bulk, variance modelling, correlation, marker
  detection) on the 10x 1.3M brain-cell dataset (`EH1039`) across three backends
  (in-memory `dgCMatrix`, `HDF5Array`, `DuckDBMatrix`), and writes
  `benchmark_results.rds`.
- **`make_timings_table.R`** — helpers the vignette `source()`s to render the
  results (`load_vignette_timings`, `make_timings_table`, `timings_config_note`).
  If `benchmark_results.rds` is absent, the vignette shows a short note instead.

## Regenerating the results

```sh
# full run (12,500 cells from EH1039) — needs ExperimentHub, scran, scuttle, HDF5Array
Rscript inst/scripts/run_scran_scuttle_benchmarks.R
# then copy benchmark_results.rds into inst/scripts/ and commit it
```

Environment variables (all optional):

| Variable | Default | Meaning |
|----------|---------|---------|
| `BENCH_NCELLS` | 12500 | cell subset from EH1039 |
| `BENCH_CORES` | available − 1 | DuckDB thread budget |
| `BENCH_SYNTHETIC` | unset | use a small synthetic matrix instead of EH1039 |

To smoke-test the harness without the ~2 GB download:

```sh
BENCH_SYNTHETIC=1 BENCH_NCELLS=400 Rscript inst/scripts/run_scran_scuttle_benchmarks.R
```

In-memory `dgCMatrix` and `HDF5Array` run single-threaded; `DuckDBMatrix` autotunes
DuckDB's internal threads up to the core budget. The per-run configuration is
recorded in `attr(results, "config")` and shown beneath the vignette table.
