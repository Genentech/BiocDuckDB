# BiocDuckDB

*High-level Bioconductor integration for DuckDB-backed data — read and write Experiment objects as Parquet, analyze them on disk.*

## Overview

`BiocDuckDB` is the integration layer of the **BiocDuckDB** suite. It brings
together [DuckDBDataFrame](https://bioconductor.org/packages/DuckDBDataFrame),
[DuckDBArray](https://bioconductor.org/packages/DuckDBArray), and
[DuckDBGRanges](https://bioconductor.org/packages/DuckDBGRanges) to provide:

- **`writeParquet()` / `readParquet()`** — serialize a `SummarizedExperiment`,
  `SingleCellExperiment`, `MultiAssayExperiment`, or `MultiAssaySpatialExperiment`
  to a self-describing directory of columnar Parquet, and read it back as a
  **DuckDB-backed** experiment with the same Bioconductor API (assays become
  `DuckDBArray`, `rowData`/`colData` become `DuckDBDataFrame`, `rowRanges` becomes
  `DuckDBGRanges`).
- **SQL-optimized [scuttle](https://bioconductor.org/packages/scuttle) and
  [scran](https://bioconductor.org/packages/scran)** — QC metrics, normalization,
  pseudo-bulk aggregation, variance modelling, and marker detection implemented as
  DuckDB queries, so they run directly on a `DuckDBMatrix` without realizing it.

The on-disk format is a [Frictionless Data Package](https://frictionlessdata.io/),
readable by other tools and languages, not just R.

## Installation

```r
# once available from Bioconductor:
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("BiocDuckDB")
```

## Quick start

```r
library(BiocDuckDB)
library(SummarizedExperiment)

data(airway, package = "airway")

# write to Parquet, read back as DuckDB-backed
path <- file.path(tempdir(), "airway_se")
writeParquet(airway, path)
se <- readParquet(path)

assay(se, "counts")                 # a DuckDBArray, on disk
colSums(assay(se[1:1000, ], "counts"))   # pushed to SQL
```

Many `scran`/`scuttle` steps run directly on the DuckDB-backed matrix:

```r
library(scuttle)
perCellQCMetrics(assay(se, "counts"))    # SQL-optimized, no realization
```

## The filter, realize, analyze pattern

For data larger than memory: filter and summarize on disk with DuckDB, realize the
manageable subset into memory, run the standard pipeline, and persist results back
to Parquet. See the *Introduction to BiocDuckDB* vignette.

## Documentation

- **Introduction to BiocDuckDB** — the `writeParquet`/`readParquet` workflow,
  storage layout, operations pushed to SQL, genomic coordinates, single-cell data,
  and the filter-realize-analyze pattern (`vignettes/BiocDuckDB.Rmd`).
- **Benchmarking BiocDuckDB** — the SQL-optimized `scran`/`scuttle` methods on a
  `DuckDBMatrix` vs in-memory and `HDF5Array`, from precomputed full-scale results
  (`vignettes/BiocDuckDB-benchmark.Rmd`).

## License

MIT License. Copyright Genentech, Inc., 2026. See `inst/COPYRIGHTS` for bundled
third-party schema components.
