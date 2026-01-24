# BiocDuckDB

## High-Level Bioconductor Integration for DuckDB-Backed Data

BiocDuckDB is the integration layer that brings together `DuckDBDataFrame`, `DuckDBArray`, and `DuckDBGRanges` to provide seamless Parquet I/O and SQL-optimized implementations of single-cell analysis methods from `scran` and `scuttle`.

### Why BiocDuckDB?

BiocDuckDB provides the missing pieces for end-to-end Bioconductor workflows:

- **`writeParquet()` / `readParquet()`**: Serialize complete *Experiment objects with full metadata preservation
- **SQL-optimized scuttle**: QC metrics, normalization, and aggregation run directly on DuckDB
- **SQL-optimized scran**: Variance modeling, marker detection, and correlation analysis in SQL
- **No code changes**: Drop-in replacement for in-memory workflows

## Performance Highlights

### Single-Cell Analysis

On a 12,500 cell subset of the 10x Genomics 1.3M Brain Cell Dataset, BiocDuckDB significantly outperforms HDF5Array:

| Operation | DuckDB (sec) | HDF5Array (sec) | Speedup |
|-----------|--------------|-----------------|---------|
| **scuttle integration** |
| `perCellQCMetrics` | 0.1 | 1.2 | **12x faster** |
| `perFeatureQCMetrics` | 0.1 | 2.9 | **28x faster** |
| `nexprs` | 0.1 | 1.5 | **15x faster** |
| `summarizeAssayByGroup` | 0.4 | 3.2 | **8.7x faster** |
| `normalizeCounts` | 0.1 | 1.1 | **10.6x faster** |
| **scran integration** |
| `modelGeneVar` | 0.7 | 5.4 | **7.7x faster** |
| `modelGeneVarByPoisson` | 1.2 | 4.4 | **3.7x faster** |
| `modelGeneCV2` | 0.6 | 3.8 | **5.9x faster** |
| `correlatePairs` | 0.7 | 55.2 | **79x faster** |
| `pairwiseTTests` | 3.7 | 8.7 | **2.4x faster** |
| `pairwiseBinom` | 3.5 | 9.6 | **2.7x faster** |
| `findMarkers` | 14.1 | 19.5 | **1.4x faster** |
| `scoreMarkers`* | 7.9 | 16.7 | **2.1x faster** |
| `summaryMarkerStats` | 1.4 | 11.9 | **8.8x faster** |

\* Uses normal approximation for AUC by default (~97% correlation); use `true.auc = TRUE` for exact rank-based statistics (100% correlation, ~15-25x slower).

**Notable:** `correlatePairs` is **79x faster than HDF5Array** and **57x faster than in-memory** sparse matrices due to sparse-aware SQL.

## Key Features

### Seamless Parquet I/O

```r
library(BiocDuckDB)
library(scRNAseq)

# Load example dataset
sce <- ZeiselBrainData()

# Write entire SCE to Parquet (with full metadata)
sce_path <- file.path(tempdir(), "zeisel")
writeParquet(sce, sce_path)

# Read back as DuckDB-backed SCE
sce_ddb <- readParquet(sce_path)
sce_ddb
# class: SingleCellExperiment 
# dim: 20006 3005 
# metadata(2): SuppInfo which_qc
# assays(1): counts (DuckDBMatrix)
# rownames(20006): Tspan12 Tshz1 ... mt-Rnr1 mt-Nd4l
# rowData names(1): featureType
# colnames(3005): 1772071015_C02 1772071017_G12 ... 1772063068_D01 1772066098_A12
# colData names(10): tissue group # ... level1class level2class
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(2): ERCC repeat
```

What gets written:
- ✅ All assays (as DuckDBMatrix in coordinate format)
- ✅ rowData / colData (as DuckDBDataFrame)
- ✅ rowRanges (as DuckDBGRanges, if present)
- ✅ Alternative experiments
- ✅ Metadata
- ✅ `datapackage.json` schema for reproducibility

### SQL-Optimized scuttle Methods

All `scuttle` functions work directly on `DuckDBMatrix`:

```r
library(scuttle)

# QC metrics (28x faster for features!)
qc <- perCellQCMetrics(counts(sce_ddb))
feature_qc <- perFeatureQCMetrics(counts(sce_ddb))

# Normalization (10.6x faster)
sce_ddb <- logNormCounts(sce_ddb)

# Pseudo-bulk aggregation (8.7x faster)
pseudo <- summarizeAssayByGroup(
    counts(sce_ddb),
    ids = DataFrame(
        sample = sce_ddb$sample_id,
        cluster = sce_ddb$cluster
    ),
    statistics = c("sum", "mean", "num.detected")
)
```

### SQL-Optimized scran Methods

All major `scran` functions optimized for DuckDB:

```r
library(scran)

# Variance modeling (7.7x faster)
dec <- modelGeneVar(logcounts(sce_ddb))
dec_pois <- modelGeneVarByPoisson(counts(sce_ddb))
dec_cv2 <- modelGeneCV2(counts(sce_ddb))

# Gene-gene correlation (79x faster!)
cor_res <- correlatePairs(logcounts(sce_ddb)[1:500, ])

# Differential expression
markers_t <- pairwiseTTests(logcounts(sce_ddb), 
                            groups = sce_ddb$cluster)
markers_b <- pairwiseBinom(counts(sce_ddb),
                           groups = sce_ddb$cluster)

# Marker scoring
markers <- findMarkers(logcounts(sce_ddb), 
                       groups = sce_ddb$cluster)
scores <- scoreMarkers(logcounts(sce_ddb),
                       groups = sce_ddb$cluster)
summary_stats <- summaryMarkerStats(logcounts(sce_ddb),
                                    groups = sce_ddb$cluster)
```

### MultiAssayExperiment Support

```r
# Write MAE to Parquet
mae <- MultiAssayExperiment(
    experiments = list(
        rnaseq = sce_rna,
        atacseq = se_atac
    ),
    colData = sample_metadata
)

writeParquet(mae, "multiomics")
mae_ddb <- readParquet("multiomics")

# All experiments are DuckDB-backed
experiments(mae_ddb)$rnaseq  # DuckDBMatrix assays
experiments(mae_ddb)$atacseq # DuckDBGRanges rowRanges
```

## Quick Start

### Complete SingleCellExperiment Workflow

```r
library(BiocDuckDB)
library(scRNAseq)
library(scuttle)
library(scran)
library(scater)

# Load dataset
sce <- ZeiselBrainData()

# Write to Parquet
writeParquet(sce, "zeisel_brain")

# Read as DuckDB-backed (memory-efficient)
sce_ddb <- readParquet("zeisel_brain")

# Standard workflow - all SQL-optimized!
# 1. QC
qc <- perCellQCMetrics(counts(sce_ddb))
sce_ddb <- sce_ddb[, qc$sum > 1000 & qc$detected > 500]

# 2. Normalization
sce_ddb <- logNormCounts(sce_ddb)

# 3. Feature selection
dec <- modelGeneVar(logcounts(sce_ddb))
hvgs <- getTopHVGs(dec, n = 2000)

# 4. Dimensionality reduction
# (For PCA/UMAP, materialize subset to memory)
sce_mem <- as(sce_ddb[hvgs, ], "SingleCellExperiment")
sce_mem <- runPCA(sce_mem, ncomponents = 50)
sce_mem <- runUMAP(sce_mem)

# 5. Clustering
clusters <- quickCluster(logcounts(sce_mem))

# 6. Marker detection (back to DuckDB for efficiency!)
markers <- findMarkers(logcounts(sce_ddb)[hvgs, ], 
                       groups = clusters)
```

### RangedSummarizedExperiment with DuckDBGRanges

```r
library(SummarizedExperiment)

# Create RSE with DuckDBGRanges rowRanges
rse <- SummarizedExperiment(
    assays = list(counts = ddb_counts),
    rowRanges = ddb_granges,
    colData = ddb_coldata
)

# Write to Parquet
writeParquet(rse, "rnaseq_experiment")

# Read back (fully DuckDB-backed)
rse_ddb <- readParquet("rnaseq_experiment")

# rowRanges operations (lazy SQL)
chr1_features <- rse_ddb[seqnames(rowRanges(rse_ddb)) == "chr1", ]
```

## Parquet Storage Structure

`writeParquet()` creates a structured directory:

```
experiment/
├── datapackage.json         # Frictionless Data Package metadata
├── assays/
│   └── counts/
│       └── data.parquet     # Coordinate (COO) format
├── colData/
│   └── data.parquet
├── rowData/
│   └── data.parquet
└── rowRanges/              # If present
    └── data.parquet
```

Benefits:
- **Columnar compression**: 5-10x smaller than HDF5 for sparse data
- **Cross-language**: Read by Python, Julia, Rust, etc.
- **Cloud-ready**: DuckDB reads directly from S3/GCS
- **Schema preservation**: Full metadata in `datapackage.json`

## Bioconductor Integration

BiocDuckDB provides optimized methods for:

| Package | Functions | Support |
|---------|-----------|---------|
| **MatrixGenerics** | `rowSums`, `colMeans`, `rowVars`, etc. | ✅ Via DuckDBArray |
| **scuttle** | `perCellQCMetrics`, `perFeatureQCMetrics`, `nexprs`, `summarizeAssayByGroup`, `normalizeCounts`, `logNormCounts` | ✅ SQL-optimized |
| **scran** | `modelGeneVar`, `modelGeneVarByPoisson`, `modelGeneCV2`, `correlatePairs`, `pairwiseTTests`, `pairwiseBinom`, `findMarkers`, `scoreMarkers`, `summaryMarkerStats` | ✅ SQL-optimized |
| **SummarizedExperiment** | Construction, assay access | ✅ Full support |
| **SingleCellExperiment** | Construction, assay access, altExps | ✅ Full support |
| **MultiAssayExperiment** | Construction, experiment access | ✅ Full support |

## When to Use BiocDuckDB

**Recommended for:**
- Storing/sharing SingleCellExperiment objects
- Large-scale single-cell analysis (>50k cells)
- Pseudo-bulk and grouped analyses
- Cross-language workflows (Python ↔ R)
- Cloud-based analysis pipelines
- When you want near in-memory performance without RAM constraints

**Consider alternatives when:**
- Dataset fits comfortably in memory (use `dgCMatrix`)
- You need mutable data structures (use HDF5)
- You need 10x Genomics native format (use `TENxMatrix`)
- You need array versioning (use TileDBArray)

## Architecture

BiocDuckDB builds on three foundation packages:

```
DuckDBDataFrame (tabular data, SQL operations)
    ├── DuckDBArray (DelayedArray backend)
    └── DuckDBGRanges (GenomicRanges backend)
         ↓
    BiocDuckDB (integration layer)
```

Each foundation package is independently useful:
- **DuckDBDataFrame**: Core infrastructure
- **DuckDBArray**: For matrix-centric workflows
- **DuckDBGRanges**: For genomic range workflows
- **BiocDuckDB**: For complete *Experiment workflows

## Documentation

- **[BiocDuckDB Classes](vignettes/BiocDuckDB-classes.Rmd)**: Building SummarizedExperiment, RangedSummarizedExperiment, and SingleCellExperiment with DuckDB-backed components
- **[BiocDuckDB-DelayedArray](vignettes/BiocDuckDB-DelayedArray.Rmd)**: scuttle and scran integration with comprehensive benchmarks

For foundation packages, see:
- **DuckDBDataFrame**: Tabular data and SQL operations
- **DuckDBArray**: Matrix operations and benchmarks  
- **DuckDBGRanges**: Genomic range operations

## Installation

```r
# Install all four packages
# install.packages("remotes")
remotes::install_github("your-org/DuckDBDataFrame")
remotes::install_github("your-org/DuckDBArray")
remotes::install_github("your-org/DuckDBGRanges")
remotes::install_github("your-org/BiocDuckDB")
```

## Dependencies

BiocDuckDB depends on:
- **Foundation**: DuckDBDataFrame, DuckDBArray, DuckDBGRanges
- **Bioconductor core**: BiocGenerics, S4Vectors, IRanges, GenomicRanges, S4Arrays, SparseArray, DelayedArray, MatrixGenerics
- **Bioconductor experiments**: SummarizedExperiment, SingleCellExperiment, MultiAssayExperiment
- **Single-cell**: scuttle, scran, BiocParallel, metapod
- **Data I/O**: arrow, duckdb, DBI, dplyr, dbplyr, jsonlite

## Contributing

Contributions are welcome! Please:
- Report issues through GitHub
- Include benchmarks for performance claims
- Test with real single-cell datasets
- Follow Bioconductor standards

## License

BiocDuckDB is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements

Special thanks to:
- The Bioconductor project for infrastructure and community
- Aaron Lun for scran, scuttle, and architectural guidance
- The DuckDB team for their excellent database
- The Apache Arrow project for the Parquet format
- The single-cell analysis community for feedback and testing
