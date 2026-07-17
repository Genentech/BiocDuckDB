# BiocDuckDB 0.9.26

## New features

- `writeParquet(..., cluster_by = )` threads a clustering key (`DuckDBDataFrame::zorder()` /
  `hilbert()`, or a character vector) through every `writeParquet` method to the primitive
  writers, so rows are physically ordered on write for DuckDB row-group zonemap pruning. The
  lazy `DuckDBTable`/`DuckDBDataFrame` path lowers it SQL-side (no materialization); the
  materializing `data.frame`/`DataFrame` path reorders in memory via
  `DuckDBDataFrame::clusterSort()`. Requires DuckDBDataFrame (>= 0.9.27).

# BiocDuckDB 0.9.23

## Bug fixes

- `readParquet()` now restores `factor` columns (including ordered factors) in
  the flat-table read paths, using the `categories`/`categoriesOrdered` recorded
  in the product schema on write. Previously these columns came back as
  `character`. (Factor restoration inside `GRanges`/`SelfHits` `mcols` is not yet
  covered.)

## Documentation

- Restructured the vignettes into a user-first set, replacing *Constructing
  Experiment Objects with BiocDuckDB* and *scuttle and scran integration for
  DuckDBMatrix*:
  - *Introduction to BiocDuckDB* --- the `writeParquet()`/`readParquet()` workflow,
    the Frictionless/coordinate storage layout, operations pushed to SQL, genomic
    coordinates, single-cell data, and the filter-realize-analyze pattern.
  - *Benchmarking BiocDuckDB* --- the SQL-optimized `scuttle`/`scran` methods on a
    `DuckDBMatrix` compared with in-memory and `HDF5Array`, rendered from
    precomputed full-scale results (10x 1.3M brain cells) so the vignette builds
    quickly.
- Added `inst/scripts/` with the offline benchmark generator
  (`run_scran_scuttle_benchmarks.R`) and the vignette table helpers
  (`make_timings_table.R`), following the `HDF5Array` performance-vignette
  precompute pattern.
- Rewrote the README.
