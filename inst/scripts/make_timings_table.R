## make_timings_table.R
##
## Helpers used by the "Benchmarking BiocDuckDB" vignette to render PRECOMPUTED
## scran/scuttle benchmark timings, so the vignette does not run the (multi-GB,
## multi-minute) benchmark at build time. Mirrors the HDF5Array performance
## vignette's precompute approach.
##
## Timings are produced offline by run_scran_scuttle_benchmarks.R and saved as
## inst/scripts/benchmark_results.rds. To regenerate (see that script's header):
##
##   Rscript inst/scripts/run_scran_scuttle_benchmarks.R
##   # then copy benchmark_results.rds into inst/scripts/ and commit it
##
## benchmark_results.rds is a tidy data.frame with columns
##   Operation, Backend (InMemory/HDF5Array/DuckDB), Seconds
## and attr(., "config"): a list (n_cells, n_genes, cores) describing the run.

## Load the shipped timings, or NULL if absent (so the vignette degrades cleanly).
load_vignette_timings <- function() {
    path <- system.file("scripts", "benchmark_results.rds", package = "BiocDuckDB")
    if (!nzchar(path) || !file.exists(path)) {
        return(NULL)
    }
    readRDS(path)
}

.BACKENDS <- c("InMemory", "HDF5Array", "DuckDB")

## Pivot to wide form (one row per operation), adding DuckDB-vs-HDF5Array and
## DuckDB-vs-in-memory speedup columns.
timings_wide <- function(results) {
    ops <- unique(results$Operation)
    out <- data.frame(Operation = ops, stringsAsFactors = FALSE)
    for (b in .BACKENDS) {
        out[[b]] <- vapply(ops, function(o) {
            v <- results$Seconds[results$Operation == o & results$Backend == b]
            if (length(v) == 0L) NA_real_ else v[[1L]]
        }, numeric(1))
    }
    out[["vs_HDF5"]] <- round(out[["HDF5Array"]] / out[["DuckDB"]], 1)
    out[["vs_InMemory"]] <- round(out[["InMemory"]] / out[["DuckDB"]], 1)
    out
}

## Render the timings as a knitr::kable, or a short note if absent.
make_timings_table <- function(caption = NULL, results = load_vignette_timings()) {
    if (is.null(results)) {
        cat("_Precomputed benchmark results are not bundled in this build. ",
            "Generate them with `inst/scripts/run_scran_scuttle_benchmarks.R` ",
            "(see that script's header)._\n", sep = "")
        return(invisible(NULL))
    }
    knitr::kable(timings_wide(results), digits = 2, caption = caption,
                 col.names = c("Operation", "In-memory (s)", "HDF5Array (s)",
                               "DuckDB (s)", "vs HDF5Array (x)", "vs in-memory (x)"))
}

## Emit the run configuration as a Markdown paragraph (for reproducibility).
timings_config_note <- function(results = load_vignette_timings()) {
    cfg <- if (is.null(results)) NULL else attr(results, "config")
    if (is.null(cfg)) {
        return(invisible(NULL))
    }
    cat(sprintf(paste0(
        "Configuration: %s genes x %s cells, %d-core budget. ",
        "**In-memory:** `dgCMatrix`. **HDF5Array:** 10x/HDF5 backend. ",
        "**DuckDB:** `DuckDBMatrix` over Parquet (autotuned threads). All backends ",
        "run the same scran/scuttle generics.\n"),
        format(cfg$n_genes, big.mark = ","),
        format(cfg$n_cells, big.mark = ","), cfg$cores))
}
