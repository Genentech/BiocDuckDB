# BiocDuckDB 0.99.20

## Bug fixes

- `.compute_true_auc_sql_DuckDBMatrix()` (the `scoreMarkers(..., true.auc = TRUE)`
  path) interpolated group labels, arbitrary user-supplied text from the
  `groups` factor's levels, straight into raw SQL (`WHERE group_label IN
  ('...', '...')`) with no escaping. A label containing a single quote (not an
  exotic case for a cell-type or condition name) produced malformed SQL.
  Found during an independent code-quality review of the package's `for`
  loops; every other raw-SQL-building spot in this file only interpolates
  internal identifiers or integers, never free-text user data, so this was the
  one place that departed from that safer pattern. Fixed by quoting both
  group labels with `DBI::dbQuoteString()` before interpolation. Verified
  with a group label containing an apostrophe.

## Documentation

- Fixed `` `r CRANpkg("airway")` `` to `` `r Biocpkg("airway")` `` in the
  introduction vignette; `airway` is a Bioconductor experiment-data package,
  not a CRAN package, so the old macro rendered a broken link.
- Fixed two comma splices in the introduction and benchmarking vignettes.
- Added a one-sentence explanation of `DuckDBMatrix()`'s `datacol`/`keycols`
  arguments where the constructor first appears in the introduction vignette,
  since neither is explained anywhere else in the vignette.

## Changes

- Vectorized the Monte-Carlo simulation loop in
  `.generate_poisson_values_DuckDBMatrix()` (used by `modelGeneVarByPoisson`'s
  trend fitting): one batched `rpois()`/`rnbinom()` draw across all simulation
  points instead of up to `npts` (default 1000) separate calls.
- Removed a dead accumulator variable (`gene_means_by_block`) in
  `.compute_blocked_stats_DuckDBMatrix()` that was written but never read.
- Collapsed 4 per-group/per-column `for` loops in `numDetectedAcrossFeatures()`,
  `sumCountsAcrossFeatures()`, and `summarizeAssayByGroup()` (sum and detected
  branches) in `DuckDBMatrix-scuttle.R` into one-line matrix broadcasts,
  matching the `t(t(...))` idiom already used elsewhere in the same file.

# BiocDuckDB 0.99.19

## Bug fixes

- The categorical value-label decoder in `.readParquetResource()` (used to
  reconstruct dictionary-encoded key columns, e.g. `sample_map`'s `assay`
  column) is built with `sapply()` specifically because its output must
  retain the names `FUN` returns (`value -> label`). During an independent
  code-quality review, this call was one of several converted to `vapply()`
  to address BiocCheck's `sapply()` NOTE; `vapply()` never takes names from
  what `FUN` returns when `FUN.VALUE` is an unnamed scalar template, so the
  conversion silently produced an unnamed vector. This corrupted keycol
  decoding, which surfaced as `sampleMap` reconstruction failures
  (`"not all ExperimentList samples are found in the sampleMap"`) on
  `writeParquet()`/`readParquet()` round-trips of spatial
  `MultiAssaySpatialExperiment` objects. Caught by the full test suite before
  being committed; reverted to `sapply()` with a comment explaining why, and
  left as the one deliberate exception BiocCheck still flags.

## Documentation

- Fixed roughly 69 `\link{}` cross-references across 16 man pages
  (`readParquet`, `writeParquet`, `DuckDBDualSubset`, the
  `DuckDBMatrix`/`DuckDBDataFrame`/`SingleCellExperiment` method families, and
  the `MultiAssaySpatialExperiment` query helpers) that pointed at classes or
  functions defined in other packages (`DuckDBArray`, `DuckDBDataFrame`,
  `DuckDBGRanges`, `MultiAssaySpatialExperiment`, `S4Vectors`) without a
  package anchor, which `R CMD check --as-cran` flags as a NOTE.

## Changes

- Added `S4Arrays` to `Suggests:`; it was already an undeclared dependency of
  two test files (`test-writeParquet-coord-append.R`,
  `test-writeParquet-roundtrip-append.R`).
- Removed 6 dead internal functions with no remaining call sites:
  `.extractIntType`, `.fieldToDuckDBCast` (`fieldtypes.R`);
  `.filterSpatialLayerByInstances`, `.filterSpatialLayersByInstances`,
  `.spatialExtent` (`MultiAssaySpatialExperiments-internals.R`); `.prefixSeq`
  (`writeParquet.R`).
- Collapsed 4 `paste()`/`paste0()`-built strings passed to `stop()`/`warning()`
  into single string literals (`DuckDBMatrix-scran.R`,
  `writeStreamingResource.R`), matching BiocCheck's preference for a single
  literal over a dynamically assembled condition message.
- Converted 8 of 9 flagged `sapply()` calls in `readParquet.R` to `vapply()`
  for a guaranteed return type and length; the 9th is documented above.

# BiocDuckDB 0.99.18

## Bug fixes

- `writeParquet.GenomicRanges()` silently wrote positional numbers
  ("1", "2", ...) instead of a named `GenomicRanges`' real feature names as
  the `__name__` display column, whenever the installed `GenomicRanges` was
  1.65.2 or later. `GenomicRanges` 1.65.2 rewrote `as.data.frame.GenomicRanges()`:
  it used to default `row.names` to `names(x)` when not supplied; the new
  version hardcodes `row.names = NULL` internally (an explicitly-passed
  `row.names` argument is silently discarded too) and puts `names(x)` in a
  new `names` column instead. `writeParquet.GenomicRanges()`'s
  `df <- as.data.frame(x, optional = TRUE)` relied on the old convention to
  carry real names into `df`'s row names, which then became the written
  `__name__` column; with the new convention, `df`'s row names were the
  default numeric sequence instead. Reading a `RangedSummarizedExperiment`
  written this way surfaced as `rowRanges()` names coming back wrong, and
  then, one step further, name-based row indexing (`x[names(x)]`) failing
  outright with a DuckDB `Conversion Error` trying to match the real names
  against an integer keycol. Root-caused by installing the
  exact upstream `GenomicRanges` commit and reproducing the original failure
  byte-for-byte; fixed by explicitly setting `rownames(df) <- names(x)`
  (and dropping the new `names` column) after the conversion, which works
  against both old and new `GenomicRanges`. `writeParquet.GenomicRangesList()`
  was checked and already set `rownames()` explicitly, so it was unaffected.
- `checkDuckDBGRanges()` (this package's own test helper, duplicated from
  `DuckDBGRanges`) had the same reliance on the old `as.data.frame.GenomicRanges()`
  convention in its own reference comparison; fixed the same way. See
  `DuckDBGRanges` 0.99.7's `NEWS.md` for the sibling fix to its copy of this
  helper.

# BiocDuckDB 0.99.17

## Bug fixes

- `.readParquetMAE()` silently misbehaved on two edge cases in a
  `multi_assay_experiment` package: a missing `subjects` resource collapsed
  `file.path(path, NULL)` to `character(0)`, reaching `DuckDBDataFrame()`
  with an opaque downstream error instead of a clear one; a missing
  `sample_map` resource indexed into a `NULL` resource the same way.
  `subjects` now errors immediately with a message naming the missing
  resource. A missing `sample_map` is now handled by deriving a real map,
  the same way `MultiAssayExperiment()`'s own constructor does when its
  `sampleMap` argument is omitted (matching column names against `colData`
  row names) rather than falling back to an empty `DataFrame`, since an
  empty `sampleMap` only satisfies `MultiAssayExperiment`'s own
  `.checkSampleNames()` validity check when every experiment is also empty.
  That derivation now lives in `.newMAE()` (this package's own internal
  constructor, called from exactly one place) rather than in
  `.readParquetMAE()` itself, so any future caller passing `sampleMap = NULL`
  gets the same behavior for free, without `.newMAE()` taking on
  `MultiAssayExperiment()`'s own `.harmonize()`/full `validObject()` cost,
  which this package's writer output does not need paying again on read.

# BiocDuckDB 0.99.16

## Bug fixes

- `writeStreamingResource()` created its `_INCOMPLETE` marker *after* part 0
  had been written, so a crash during part 0, the very case the marker
  exists to catch. It left no marker at all, contradicting the documented
  behavior. It is now created immediately before part 0, along with the
  directory to hold it.

# BiocDuckDB 0.99.15

- Reading a `MultiAssaySpatialExperiment` loaded every raster-backed label
  twice: correctly into `spatialLabels()`, and again into `spatialImages()`
  under the unstripped name `sample_labels_<nm>`. Only an array-valued label
  takes the `spatial_label_coord` branch on write; a raster-backed one is
  written with the same `spatial_raster_ref` layout an image uses, so
  selecting the image resources by layout alone could not tell the two apart.
  The label branch already narrowed its own selection by name; the image
  branch now does the same.

- `readParquet()` could not read a package declaring
  `model = "experiment_list"`; it failed with `subscript out of bounds`. The
  dispatch passed the whole `datapackage.json` manifest to
  `.readParquetExps()`, which takes the *resource list* rather than the
  manifest (its two other callers, the `altExps` path and
  `.readParquetMAE()`, both pass a filtered resource list). It therefore
  iterated the manifest's top-level entries and tried to read a `name` field
  off `$schema`. It is now handed `package[["resources"]]`.

  `experiment_list` is one of the reader values named in the profile's
  `model` description, and `writeParquet()` has had an `ExperimentList`
  method throughout. The gap went unnoticed because the round-trip suite
  explicitly skipped `ExperimentList` (a note to that effect sat in
  `test-parquet-roundtrip.R`); it is now covered.

# BiocDuckDB 0.99.14

## Testing

- Loading the `airway` Suggests package in `tests/testthat/setup.R` is now
  conditional, and the airway-dependent test files
  (`test-DuckDBMatrix-scuttle.R`, `test-DuckDBMatrix-scran.R`) and test block
  (`test-parquet-roundtrip.R`) call `skip_if_not_installed("airway")`.
  Previously the unconditional `data(airway, ...)` call in `setup.R` would
  fail and abort the entire test suite (not just the airway-dependent tests)
  when `airway` was not installed. Same class of issue as found in the
  DuckDBArray Bioconductor review.
- `test-MultiAssaySpatialExperiment.R` calls `skip_if_not_installed("sf")`,
  the Suggests package its fixtures build geometries with.

# BiocDuckDB 0.99.13

## Follow-up review changes

- Updated the sole caller of the shared internal generic (a round-trip test) to
  the undotted name after the Bioconductor-review rename in DuckDBDataFrame
  (`.has_row_number` -> `has_row_number`; guideline: functions starting with `.`
  should not be exported). Requires the renamed version of DuckDBDataFrame.

## Documentation

- Replaced em dashes with commas or colons in the vignettes.

# BiocDuckDB 0.99.12

## Documentation

- Added `URL` and `BugReports` fields to DESCRIPTION.
- Added a package-level man page (`?BiocDuckDB`).
- Removed the redundant `library(BiocStyle)` call from the vignettes (the
  `BiocStyle::html_document` output and `::`-qualified helpers apply the style
  without loading the package). Applied suite-wide for consistency with the
  Bioconductor review of DuckDBDataFrame.

# BiocDuckDB 0.99.11

## Enhancements

- `readParquet()` now assembles `multi_assay_experiment` and
  `multi_assay_spatial_experiment` objects directly with
  `S4Vectors::new2(check = FALSE)` instead of calling the
  `MultiAssayExperiment()` / `MultiAssaySpatialExperiment()` constructors. Those
  constructors re-derive sample-name set equality at construction time via a
  locale-collation sort over the full, cell-grain `sampleMap`, which dominates
  read time on large products (tens of millions of `sampleMap` rows -- e.g. a
  cell-resolved spatial atlas took ~22 min to read, over 60% of it in that
  redundant sort). The datapackage already carries `writeParquet()`'s write-time
  integrity guarantee (a DuckDB anti-join proving every `sampleMap` `primary` is
  a subjects rowname), so the constructor-time re-check is redundant; the
  assembled object still passes `validObject()`. The direct assembly also skips
  harmonization, so unreferenced `colData` rows are retained rather than dropped
  on read (write with `subset_subjects_to_referenced` upstream to avoid them).

## Testing

- The `multi_assay_experiment` and `multi_assay_spatial_experiment` round-trip
  tests now assert `validObject()` on every object returned by `readParquet()`.
  Because the read path assembles with `new2(check = FALSE)` (no constructor-time
  validity), these explicit checks guard against a `readParquet()`-assembled
  object being structurally invalid.

# BiocDuckDB 0.99.10

## Enhancements

- `writeParquet()` now serializes a `MultiAssaySpatialExperiment`'s `spatialMap`
  as a monomorphic bridge over a conformed **element-instance registry**
  (`spatial_element_registry`), so the observation-to-spatial-element association
  becomes real, single-target foreign keys instead of an opaque polymorphic
  reference. A new `spatial_element_registry` resource enumerates every
  points/shapes instance across all layers with an integer `__element__` spine;
  each typed spatial layer foreign-keys that spine; and `spatial_map` declares
  two foreign keys -- the element side (`__element__` to `spatial_element_registry`)
  and the observation side (`(assay, colname)` to `sample_map`, which encodes
  the spatialMap validity rule as a real composite FK). The `__element__` spine
  is a serialization-internal column: `readParquet()` reconstructs the
  object-model `spatialMap` and the spatial layers without it, so the round-trip
  is unchanged.

## Bug fixes

- `writeParquet()` emitted a `coord_array` index column's category `label` as a
  one-element array rather than a string, so a datapackage carrying an assay with
  dimnames failed validation against the bundled Frictionless profile. The label
  is now a string, and such packages conform.

# BiocDuckDB 0.99.9

## Enhancements

- `writeParquet()` now declares dimension foreign keys for neighbor graphs and
  reducedDim embeddings, so a `datapackage.json` losslessly serializes the object
  model's relational structure rather than only the assay's. A `colPairs` /
  `rowPairs` `graph_edges` resource declares both edge endpoints (`from` / `to`)
  as foreign keys into the sample / feature dimension's crossing key (a
  role-playing self-reference), and a `reducedDims` / `rowLoadings`
  `embedding_table` resource declares its key column as a foreign key into its
  dimension table (an outrigger). The additions are additive and reader-neutral:
  the reader reconstructs graphs and embeddings from `graphEdges` / the Parquet
  schema as before, and the emitted package still validates against the bundled
  Frictionless profile.

# BiocDuckDB 0.99.8

## New features

- `writeStreamingResource()` now guards a multi-part write with an in-progress
  marker: an `_INCOMPLETE` file is dropped into the resource directory once part
  0 creates it and removed once the stream finishes. A completed resource is
  therefore a clean directory of only Parquet parts that reads normally, while an
  interrupted write leaves the marker behind; because the marker is not a Parquet
  file, reading the incomplete directory fails loudly instead of silently
  returning the partial parts. This makes "readable only when complete" hold by
  construction with no change to the reader.

## Enhancements

- The flat `data.frame` Parquet write path now pins the row-group size
  (`arrow::write_parquet(chunk_size = ...)`, default 491520 rows, matching the
  coordinate-array writer) instead of relying on the Arrow default. A large or
  streamed flat resource is written with bounded row groups, so range predicates
  prune at row-group granularity and reads parallelize across groups; a caller
  may still override `chunk_size`.

# BiocDuckDB 0.99.7

## New features

- `writeStreamingResource()` streams a larger-than-memory table to a single
  flat, multi-part Parquet resource by pulling it block-by-block from a producer
  callback (`blocks(i)` returning the i-th `data.frame` or `NULL` when
  exhausted). It owns the streaming bookkeeping that large producers otherwise
  hand-roll on top of `writeParquet()`: the running row `offset`, the `part`
  index and zero-padded `part_digits`, the `append` flag, the positional
  `dimtbl` slice per block, and `index_max` threading so the `__index__` column
  stays one consistent (possibly > 2^31) integer type across parts. It also
  performs the post-write integrity checks — a fatal partition-alignment guard
  when the streamed row count disagrees with the dimension table, a coverage
  warning against `expected_rows`, and a narrowing-floor warning when a small
  part 0 is followed by more parts without `index_max`. The returned descriptor
  (with `n_rows`) is ready for `writeDatapackage()`.

# BiocDuckDB 0.99.6

## Bug fixes

- Reading a `graph_edges` resource no longer materializes the `__index__` keycol.
  `.readParquetGraphEdges()` now constructs the graph with a `row_number` key
  (`keycol = NULL`) instead of the schema's `__index__` keycol. A level-discovered
  keycol pulls every distinct edge id into R and sorts it, which overflowed
  `bit64`'s radixsort ("long vectors not supported") for a graph with more than
  ~2.1e9 edges. A graph edge is identified by its `from`/`to` endpoints, not its
  row position, so the row-number key reconstructs the same graph without the
  scan.

# BiocDuckDB 0.99.5

- `.writeDataFrameParquet()` now fails loudly when a part's `__index__` would
  exceed the 32-bit range but `index_max` was not supplied, instead of silently
  writing a floating-point (`float64`) key spine. Producers streaming past 2^31
  rows must declare `index_max` (or `Inf`) so the column is typed `int64`
  consistently across parts.

- `writeDatapackage()` validates its inputs at the write seam (every resource
  must be a descriptor with a single-string `name` and `path`; names must be
  unique), so a malformed manifest fails fast at write time rather than
  surfacing as an obscure error at read time.

# BiocDuckDB 0.99.4

## Bug fixes

- `writeParquet()` (flat / data-frame path) can now write a resource with more
  than ~2.1e9 rows without the `__index__` column overflowing 32-bit integers. A
  new optional `index_max` argument declares the resource's total index range;
  when it exceeds the 32-bit limit the `__index__` column is typed as a 64-bit
  integer up front (part 0 included) so every streamed append part shares one
  type, instead of narrowing part 0 to `int32` and overflowing later parts on
  cast (mirrors the coord-array `max_dim` typing). Pass `index_max = Inf` when
  the total row count is unknown before streaming (e.g. graph edges). Small
  resources are unaffected --- the index still narrows to the smallest integer
  type. Requires the companion `DuckDBDataFrame` `validateAppendOffset()` fix.

# BiocDuckDB 0.99.3

## New features

- Added the exported `writeDatapackage()` function, which assembles and writes a
  Frictionless `datapackage.json` envelope from a list of resource descriptors.
  The experiment-level `writeParquet()` methods (`SummarizedExperiment`,
  `MultiAssayExperiment`) now single-source their manifest assembly through it,
  and producers that build resources incrementally --- streaming a dataset too
  large to hold in memory, or promoting from another store --- can emit a
  conformant manifest without reconstructing an in-memory Bioconductor object.
  `NULL` descriptors (returned by append/streaming parts) are dropped, so
  accumulated resource lists can be passed straight through. This is the write
  half of the ingest contract; the read half is `readParquet()` together with
  the `DuckDBMatrix()`/`DuckDBArray()`/`DuckDBTable()` constructors, which attach
  existing Parquet in place.

## Documentation

- Documented the storage layout as a targetable Frictionless contract in the
  package vignette (the new "Targeting the storage contract" section), covering
  `writeDatapackage()` for assembling a manifest and the DuckDB-backed
  constructors for attaching existing coord-array Parquet in place.

# BiocDuckDB 0.99.2

## Documentation

- Added `\value` sections (roxygen `@return`) to the `DuckDBDataFrame-spatial`,
  `DuckDBDualSubset-class`, `DuckDBMatrix-scran`, `DuckDBMatrix-scuttle`, and
  `MultiAssaySpatialExperiment-spatial` man pages, documenting the values
  returned by the spatial query, dual-subset, scran, scuttle, and MASE spatial
  I/O methods. Resolves the `R CMD BiocCheck` "missing \value" WARNING.

# BiocDuckDB 0.99.1

## Bug fixes

- The spatial-map query functions (`linkSpatialMap()`, `spatialViews()`,
  `spatialCoordinateSystems()`, `validateSpatialMap()`) now have runnable
  examples backed by a minimal spatial `MultiAssaySpatialExperiment` fixture
  bundled under `inst/extdata/spatial_mase/` (generated by
  `inst/scripts/make-spatial-mase-fixture.R`). This raises the share of
  exported-object man pages with runnable examples above the 80% BiocCheck
  threshold, resolving the `R CMD BiocCheck` ERROR.

# BiocDuckDB 0.9.29

## Bug fixes

- `writeParquet()` now right-sizes the internal integer index columns of
  data-frame resources (the `__index__`/`__sample__`/`__feature__` key spine and
  graph-edge `from`/`to`) to their smallest Arrow integer type. This narrowing
  was previously unreachable: it was gated on `!flat_part`, but
  `setupFlatParquetWrite()` always resolves `flat_part = TRUE`, so data-frame
  index columns were written as `int32`. The coord-array assay path already
  narrowed, so only the tabular resources were affected.
  - Scoped to internal index columns to match `scibis`
    (`_optimize_integer_columns` / `_narrow_index`): user data columns and
    Hive partition-group columns keep their declared dtype, so a write→read
    round-trip is dtype-preserving.
  - On append, each part's index columns are pinned to the existing part-0
    schema (via `reconcileParquetSchema()`, whose result was previously
    discarded) so every part of a streamed resource shares one integer width.
  - The recorded `arrowType` in each resource's schema now reflects the narrowed
    width, matching the physical Parquet.

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
