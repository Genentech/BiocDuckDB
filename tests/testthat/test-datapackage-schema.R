# Conformance: datapackage.json emitted by writeParquet() must validate against
# the bundled BiocDuckDB Frictionless profile (inst/schema/biocduckdb-profile.json),
# which extends the unmodified Frictionless Data Package v2.0 base spec. AJV
# resolves the relative $ref to the sibling datapackage.json offline.
# Run: library(BiocDuckDB); library(testthat); source("setup.R"); source("test-datapackage-schema.R")

library(SingleCellExperiment)
library(MultiAssayExperiment)

.profileValidator <- function() {
    skip_if_not_installed("jsonvalidate")
    profile <- system.file("schema", "biocduckdb-profile.json", package = "BiocDuckDB")
    if (!nzchar(profile) || !file.exists(profile)) {
        skip("bundled schema profile not found (package not installed with inst/schema)")
    }
    jsonvalidate::json_schema$new(profile, engine = "ajv")
}

.validateDataPackage <- function(validator, dir) {
    dp <- file.path(dir, "datapackage.json")
    expect_true(file.exists(dp))
    json <- paste(readLines(dp, warn = FALSE), collapse = "\n")
    validator$validate(json)
}

test_that("SummarizedExperiment datapackage.json conforms to the profile", {
    validator <- .profileValidator()

    set.seed(200)
    ngenes <- 40L
    ncells <- 20L
    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = DataFrame(gene_type = sample(c("protein_coding", "lncRNA"),
                                               ngenes, replace = TRUE)),
        colData = DataFrame(batch = sample(c("A", "B"), ncells, replace = TRUE))
    )

    tmpdir <- tempfile()
    writeParquet(se, tmpdir)
    expect_true(.validateDataPackage(validator, tmpdir))
    unlink(tmpdir, recursive = TRUE)
})

test_that("SingleCellExperiment with embeddings conforms (fixed-length array field)", {
    validator <- .profileValidator()

    set.seed(201)
    ngenes <- 40L
    ncells <- 20L
    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    sce <- SingleCellExperiment(assays = list(counts = counts))
    reducedDim(sce, "PCA") <- matrix(rnorm(ncells * 5L), nrow = ncells, ncol = 5L)
    reducedDim(sce, "UMAP") <- matrix(rnorm(ncells * 2L), nrow = ncells, ncol = 2L)

    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)

    # The embedding tables carry array-typed fields with equal min/max length;
    # this exercises the only field-level extension (arrayItem) under the profile.
    expect_true(.validateDataPackage(validator, tmpdir))
    unlink(tmpdir, recursive = TRUE)
})

test_that("MultiAssayExperiment datapackage.json conforms to the profile", {
    validator <- .profileValidator()

    set.seed(202)
    ngenes <- 30L
    ncells <- 15L
    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))
    se <- SummarizedExperiment(assays = list(counts = counts))

    mae <- MultiAssayExperiment(experiments = ExperimentList(rna = se))

    tmpdir <- tempfile()
    writeParquet(mae, tmpdir)
    expect_true(.validateDataPackage(validator, tmpdir))
    unlink(tmpdir, recursive = TRUE)
})

test_that("MultiAssayExperiment sample_map declares the primary -> subjects foreign key (ADR-053)", {
    set.seed(203)
    ngenes <- 10L
    ncells <- 6L
    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))
    se <- SummarizedExperiment(assays = list(counts = counts))
    mae <- MultiAssayExperiment(experiments = ExperimentList(rna = se))

    tmpdir <- tempfile()
    writeParquet(mae, tmpdir)
    dp <- jsonlite::fromJSON(file.path(tmpdir, "datapackage.json"),
                             simplifyVector = FALSE)
    get_res <- function(nm) Find(function(r) identical(r[["name"]], nm),
                                 dp[["resources"]])
    sm <- get_res("sample_map")
    subj <- get_res("subjects")
    expect_false(is.null(sm))
    expect_false(is.null(subj))

    fks <- sm[["schema"]][["foreignKeys"]]
    expect_equal(length(fks), 1L)
    expect_identical(fks[[1L]][["fields"]], "primary")
    expect_identical(fks[[1L]][["reference"]][["resource"]], "subjects")
    # the FK targets whatever key column 'subjects' actually emitted
    expect_identical(fks[[1L]][["reference"]][["fields"]],
                     subj[["schema"]][["primaryKey"]])

    unlink(tmpdir, recursive = TRUE)
})
