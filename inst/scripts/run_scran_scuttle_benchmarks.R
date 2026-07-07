#!/usr/bin/env Rscript
# Benchmark for the "Benchmarking BiocDuckDB" vignette.
#
# Times the SQL-optimized scuttle/scran methods BiocDuckDB implements for
# DuckDBMatrix against the same generics on an in-memory dgCMatrix and on
# HDF5Array, on the 10x Genomics 1.3M brain-cell dataset (EH1039):
#
#   perCellQCMetrics, perFeatureQCMetrics   -- scuttle QC
#   summarizeAssayByGroup                   -- pseudo-bulk (GROUP BY)
#   normalizeCounts                         -- scuttle normalization
#   modelGeneVar                            -- scran variance modelling
#   correlatePairs                          -- scran pairwise correlations (top HVGs)
#   pairwiseTTests, findMarkers             -- scran marker detection
#
# In-memory dgCMatrix and HDF5Array are single-threaded here; DuckDBMatrix
# autotunes DuckDB's internal threads. Each timing is wrapped so an
# unsupported/failed op records NA rather than aborting the run.
#
# Writes benchmark_results.rds: a tidy data.frame (Operation, Backend, Seconds)
# with attr(., "config") (n_cells, n_genes, cores). The vignette renders it via
# inst/scripts/make_timings_table.R.
#
# Usage:  Rscript run_scran_scuttle_benchmarks.R
# Env:    BENCH_NCELLS     cell subset (default 12500)
#         BENCH_CORES      DuckDB thread budget (default = available - 1)
#         BENCH_SYNTHETIC  if set, use a small synthetic matrix instead of EH1039
#                          (smoke-test the harness without the ~2 GB download)

suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
    library(BiocDuckDB)
    library(DuckDBArray)
    library(Matrix)
    library(scuttle)
    library(scran)
    library(BiocParallel)
    library(DBI)
})

## ---- configuration -----------------------------------------------------------
Sys.unsetenv("IS_BIOC_BUILD_MACHINE")   # offline benchmark; honor the real budget

available_cores <- function() {
    for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")) {
        x <- suppressWarnings(as.integer(Sys.getenv(v, "")))
        if (length(x) == 1L && !is.na(x) && x > 0L) return(x)
    }
    np <- tryCatch(as.integer(system("nproc", intern = TRUE, ignore.stderr = TRUE)),
                   error = function(e) NA_integer_, warning = function(w) NA_integer_)
    if (length(np) == 1L && !is.na(np) && np > 0L) return(np)
    max(1L, parallel::detectCores())
}
n_cells <- as.integer(Sys.getenv("BENCH_NCELLS", "12500"))
cores <- as.integer(Sys.getenv("BENCH_CORES",
                                as.character(max(1L, available_cores() - 1L))))
BPPARAM <- SerialParam()
setAutoBlockSize(2^30)
setAutoBlockShape("scale")
setAutoBPPARAM(BPPARAM)

## ---- data --------------------------------------------------------------------
cat("Preparing backends...\n")
if (nzchar(Sys.getenv("BENCH_SYNTHETIC"))) {
    set.seed(1L)
    ng <- 2000L
    brain_mem <- as(Matrix(rpois(ng * n_cells, 0.2), ng, n_cells, sparse = TRUE),
                    "dgCMatrix")
    rownames(brain_mem) <- paste0("Gene", seq_len(ng))
    colnames(brain_mem) <- paste0("Cell", seq_len(n_cells))
    brain_hdf5 <- writeHDF5Array(brain_mem)
} else {
    library(ExperimentHub)
    hub <- ExperimentHub()
    brain_full <- TENxMatrix(hub[["EH1039"]], group = "mm10")
    brain_hdf5 <- brain_full[, seq_len(n_cells)]
    brain_mem <- as(brain_hdf5, "dgCMatrix")
}
n_genes <- nrow(brain_mem)
cat(sprintf("  matrix: %d genes x %d cells | cores = %d\n", n_genes, n_cells, cores))

# DuckDBMatrix backend (transpose so features are columns for columnar access)
duckdb_path <- tempfile()
brain_t <- t(brain_mem)
writeParquet(brain_t, duckdb_path)
brain_ddb <- DuckDBMatrix(duckdb_path, datacol = "value",
    keycols = list(index2 = setNames(seq_len(ncol(brain_t)), colnames(brain_t)),
                   index1 = setNames(seq_len(nrow(brain_t)), rownames(brain_t))),
    dimtbls = createDimTables(brain_t))

# cap DuckDB threads at the budget
ddb_conn <- DuckDBDataFrame::acquireDuckDBConn()
try(dbExecute(ddb_conn, sprintf("SET threads = %d;", cores)), silent = TRUE)

groups <- paste0("Cluster_", sample(1:20, n_cells, replace = TRUE))

## ---- timing ------------------------------------------------------------------
elapsed <- function(expr)
    tryCatch(unname(system.time(force(expr))["elapsed"]),
             error = function(e) { message("    (failed: ", conditionMessage(e), ")"); NA_real_ })

# Log-normalized matrices are needed by the variance/marker ops.
cat("Normalizing (for variance/marker ops)...\n")
log_mem  <- normalizeCounts(brain_mem)
log_hdf5 <- normalizeCounts(brain_hdf5)
log_ddb  <- normalizeCounts(brain_ddb)
hvg <- head(order(modelGeneVar(log_ddb)$bio, decreasing = TRUE), 200)

# op -> per-backend thunk (list of three: InMemory, HDF5Array, DuckDB)
ops <- list(
    perCellQCMetrics = list(
        function() perCellQCMetrics(brain_mem),
        function() perCellQCMetrics(brain_hdf5),
        function() perCellQCMetrics(brain_ddb)),
    perFeatureQCMetrics = list(
        function() perFeatureQCMetrics(brain_mem),
        function() perFeatureQCMetrics(brain_hdf5),
        function() perFeatureQCMetrics(brain_ddb)),
    summarizeAssayByGroup = list(
        function() summarizeAssayByGroup(brain_mem, groups, statistics = "sum"),
        function() summarizeAssayByGroup(brain_hdf5, groups, statistics = "sum"),
        function() summarizeAssayByGroup(brain_ddb, groups, statistics = "sum")),
    normalizeCounts = list(
        function() normalizeCounts(brain_mem),
        function() normalizeCounts(brain_hdf5),
        function() normalizeCounts(brain_ddb)),
    modelGeneVar = list(
        function() modelGeneVar(log_mem),
        function() modelGeneVar(log_hdf5),
        function() modelGeneVar(log_ddb)),
    correlatePairs = list(
        function() correlatePairs(log_mem, subset.row = hvg),
        function() correlatePairs(log_hdf5, subset.row = hvg),
        function() correlatePairs(log_ddb, subset.row = hvg)),
    pairwiseTTests = list(
        function() pairwiseTTests(log_mem, groups = groups),
        function() pairwiseTTests(log_hdf5, groups = groups),
        function() pairwiseTTests(log_ddb, groups = groups)),
    findMarkers = list(
        function() findMarkers(log_mem, groups = groups, BPPARAM = BPPARAM),
        function() findMarkers(log_hdf5, groups = groups, BPPARAM = BPPARAM),
        function() findMarkers(log_ddb, groups = groups, BPPARAM = BPPARAM))
)

backends <- c("InMemory", "HDF5Array", "DuckDB")
rows <- list()
for (opname in names(ops)) {
    cat(sprintf("\n--- %s ---\n", opname))
    for (i in seq_along(backends)) {
        secs <- elapsed(ops[[opname]][[i]]())
        rows[[length(rows) + 1L]] <- data.frame(
            Operation = opname, Backend = backends[i], Seconds = secs,
            stringsAsFactors = FALSE)
        cat(sprintf("  %-11s %s\n", backends[i],
                    if (is.na(secs)) "NA" else sprintf("%.2f s", secs)))
    }
}
results <- do.call(rbind, rows)
attr(results, "config") <- list(n_cells = n_cells, n_genes = n_genes, cores = cores)

saveRDS(results, "benchmark_results.rds")
cat("\nSaved benchmark_results.rds\n")
print(results, row.names = FALSE)
