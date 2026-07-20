# BiocDuckDB 0.9.28

## New features

- `readParquet()` can now read a dataset from remote object storage (`s3://`,
  `gs://`, `http(s)://`, …). The `datapackage.json` sidecar is fetched over DuckDB
  `httpfs` (`jsonlite::read_json` cannot fetch object-storage schemes) after
  `DuckDBDataFrame::configureCloud()` installs/loads `httpfs` and applies any
  `s3_*` credentials; the resource reads then prune row groups over the VFS the
  same as on local disk. Local paths are unchanged.

## Changes

- `writeParquet()` refuses a remote object-storage path with an actionable error
  ("write to a local directory and upload it"). Object stores have no atomic
  directory swap and the append/refuse-to-overwrite guards rely on local file
  existence, so BiocDuckDB writes stay local + whole-object; upload the written
  directory afterwards (e.g. `aws s3 cp --recursive`).

# BiocDuckDB 0.9.27

## Changes

- The test harness pins DuckDB to a single thread
  (`options(DuckDBDataFrame.threads = 1L)` in `setup.R`) so parallel
  floating-point reductions accumulate in a fixed order and tight-tolerance
  expectations stay reproducible run-to-run. This is a test-only change; real
  sessions use all cores.

# BiocDuckDB 0.9.26

## New features

- Cross-element DuckDB query layer over a `MultiAssaySpatialExperiment` whose
  spatial layers are DuckDB-backed (`R/MultiAssaySpatialExperiment-query.R`). MASE
  stays DuckDB-free; this is BiocDuckDB's DuckDB engine operating on a MASE, using
  DuckDBSpatial (Suggests) for spatial SQL.
  - `spatialViews()` registers each spatial layer and the `spatialMap` junction as
    on-the-fly DuckDB temp views (lazy layers as views over their rendered SQL; no
    materialization).
  - `linkSpatialMap()` links assay observations to their spatial layer rows through
    the full `spatialMap` key (`assay`, `colname`, `element_type`, `region`,
    `instance_id`), returning a lazy `DuckDBDataFrame`.
  - `validateSpatialMap()` checks referential integrity via DuckDB anti-joins
    (orphan `instance_id`/`colname`, unknown layer, duplicate `instance_id`);
    reports by default, errors with `strict = TRUE`.
  - `spatialElementJoin()` runs an `ST_*` cross-element spatial join (via
    DuckDBSpatial), first aligning both elements into a common coordinate system
    through the coordinate-transform graph when `coordinate_system` is given.
  - `spatialCoordinateSystems()` and per-element transforms stored in
    `metadata(mase)$transforms` (`"<element_type>/<region>" -> {cs -> transform}`,
    RFC-5-shaped) round-trip through the standard MASE metadata annotations.
- `writeParquet(..., cluster_by = )` threads a clustering key (`DuckDBDataFrame::zorder()` /
  `hilbert()`, or a character vector) through every `writeParquet` method to the primitive
  writers, so rows are physically ordered on write for DuckDB row-group zonemap pruning. The
  lazy `DuckDBTable`/`DuckDBDataFrame` path lowers it SQL-side (no materialization); the
  materializing `data.frame`/`DataFrame` path reorders in memory via
  `DuckDBDataFrame::clusterSort()`. Requires DuckDBDataFrame (>= 0.9.27).

## Bug fixes

- Writing an in-memory geometry table (an `sf` / `DataFrame` with a geometry
  column, e.g. a `ShapesLayerList` layer) to Parquet no longer errors with "flat
  append ('append', 'part') is not supported for sf objects". Since flat writes
  became the default, every such write was a flat part 0, which the `sf` branch of
  `.writeDataFrameParquet` mistook for an append/part continuation. The guard now
  rejects only a genuine append or a subsequent part (`part > 0`), matching the
  lazy geometry write path, so a normal geometry write goes through
  `DuckDBSpatial::writeGeoParquet`.

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
