# Public ingest-contract surface: writeDatapackage() (envelope assembler) and
# addGraphMetadata() (resource-level graph-edges helper). writeDatapackage() is
# the seam the SummarizedExperiment / MultiAssayExperiment writeParquet() methods
# now single-source, and the entry point for producers that accumulate resources
# incrementally.
# Run: library(BiocDuckDB); library(testthat); source("setup.R"); source("test-writeDatapackage.R")

library(SummarizedExperiment)

.profileValidator <- function() {
    skip_if_not_installed("jsonvalidate")
    profile <- system.file("schema", "biocduckdb-profile.json", package = "BiocDuckDB")
    if (!nzchar(profile) || !file.exists(profile)) {
        skip("bundled schema profile not found (package not installed with inst/schema)")
    }
    jsonvalidate::json_schema$new(profile, engine = "ajv")
}

.readDataPackage <- function(dir) {
    jsonlite::fromJSON(file.path(dir, "datapackage.json"), simplifyVector = FALSE)
}

.featureResource <- function(name = "features") {
    list(name = name, path = name, dimension = "feature", layout = "data_frame",
         format = "parquet", mediatype = "application/vnd.apache.parquet",
         schema = list(fields = list(list(name = "__feature__", type = "string"))))
}

test_that("writeDatapackage emits a profile-conformant manifest", {
    validator <- .profileValidator()
    tmpdir <- tempfile()
    writeDatapackage("summarized_experiment", list(.featureResource()), tmpdir)

    dp <- file.path(tmpdir, "datapackage.json")
    expect_true(file.exists(dp))
    json <- paste(readLines(dp, warn = FALSE), collapse = "\n")
    expect_true(validator$validate(json))
    unlink(tmpdir, recursive = TRUE)
})

test_that("writeDatapackage sets $schema/model and drops NULL resources", {
    tmpdir <- tempfile()
    res <- list(.featureResource("features"), NULL, .featureResource("samples"))
    pkg <- writeDatapackage("summarized_experiment", res, tmpdir)

    # $schema is the leading key; model set; NULLs filtered.
    expect_identical(names(pkg)[1L], "$schema")
    expect_identical(pkg[["model"]], "summarized_experiment")
    expect_length(pkg[["resources"]], 2L)

    parsed <- .readDataPackage(tmpdir)
    expect_identical(parsed[["model"]], "summarized_experiment")
    expect_length(parsed[["resources"]], 2L)
    unlink(tmpdir, recursive = TRUE)
})

test_that("writeDatapackage includes main_exp_name / annotations only when given", {
    tmpdir <- tempfile()
    bare <- writeDatapackage("summarized_experiment", list(.featureResource()),
                             tmpdir)
    expect_null(bare[["main_exp_name"]])
    expect_null(bare[["annotations"]])

    tmpdir2 <- tempfile()
    full <- writeDatapackage("single_cell_experiment", list(.featureResource()),
                             tmpdir2, main_exp_name = "rna",
                             annotations = list(note = "hi"))
    expect_identical(full[["main_exp_name"]], "rna")
    expect_identical(full[["annotations"]], list(note = "hi"))
    unlink(c(tmpdir, tmpdir2), recursive = TRUE)
})

test_that("writeDatapackage validates its arguments", {
    expect_error(writeDatapackage(c("a", "b"), list(), tempfile()), "single")
    expect_error(writeDatapackage("m", "not-a-list", tempfile()), "list")
    expect_error(writeDatapackage("m", list(), 1L), "single")
})

test_that("writeDatapackage re-assembles a reader-valid manifest (round-trip)", {
    set.seed(11)
    counts <- matrix(rpois(30L, 5), nrow = 6L, ncol = 5L)
    rownames(counts) <- paste0("Gene", seq_len(6L))
    colnames(counts) <- paste0("Cell", seq_len(5L))
    se <- SummarizedExperiment(assays = list(counts = counts))

    dir <- tempfile()
    writeParquet(se, dir)                       # writes resources + manifest
    dp <- .readDataPackage(dir)

    # Re-emit the manifest from its own parsed resources via the public function,
    # then confirm the reader still reconstructs the object.
    writeDatapackage(model = dp[["model"]], resources = dp[["resources"]],
                     path = dir, annotations = dp[["annotations"]])
    se2 <- readParquet(dir)
    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(dim(se2), dim(se))
    expect_identical(assayNames(se2), assayNames(se))
    unlink(dir, recursive = TRUE)
})
