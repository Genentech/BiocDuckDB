# A package is read on the word of its datapackage.json. A manifest that is not
# the shape that the readers walk, or that names a place outside the package, is
# refused with a message, and a package that was written by this package reads
# without a warning.

skip_if_not_installed("SummarizedExperiment")

.safe_package <- function() {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(c(1, 0, 2, 0, 3, 4), 3, 2,
                                      dimnames = list(c("g1", "g2", "g3"), c("c1", "c2")))),
        colData = S4Vectors::DataFrame(b = 1:2, row.names = c("c1", "c2")))
    out <- tempfile()
    writeParquet(se, out)
    out
}

.edit_manifest <- function(path, FUN) {
    f <- file.path(path, "datapackage.json")
    dp <- jsonlite::read_json(f, simplifyVector = TRUE, simplifyDataFrame = FALSE,
                              simplifyMatrix = FALSE)
    jsonlite::write_json(FUN(dp), f, auto_unbox = TRUE, null = "null")
    invisible(path)
}

test_that("a package that was written reads without a warning", {
    out <- .safe_package()
    on.exit(unlink(out, recursive = TRUE))
    expect_no_warning(se <- readParquet(out))
    expect_identical(dim(se), c(3L, 2L))
})

test_that("a manifest of the wrong shape is a classed error", {
    out <- .safe_package()
    on.exit(unlink(out, recursive = TRUE))
    expect_error(readParquet(out, package = "not a manifest"),
                 class = "biocduckdb_manifest_error")
    expect_error(.validateDatapackage(list(resources = "x")), class = "biocduckdb_manifest_error")
    expect_error(.validateDatapackage(list(resources = list(list(name = "a")))),
                 "resource 1 needs a 'path'")
    expect_error(.validateDatapackage(list(resources = list(list(name = "", path = "p")))),
                 "needs a 'name'")
    expect_error(.validateDatapackage(list(resources = list(list(name = "a", path = "p"),
                                                            list(name = "b", path = "p")))),
                 "same path")
    expect_silent(.validateDatapackage(list(resources = list(list(name = "a", path = "p")))))
    expect_silent(.validateDatapackage(list(name = "no resources")))
})

test_that("a resource that names a place outside the package is refused", {
    for (bad in c("../outside/samples", "samples/../../x", "/etc/passwd", "~/x",
                  "C:/x", "\\\\server\\share", "s3://bucket/x", "a\\..\\b")) {
        expect_error(.checkResourcePath(bad), class = "biocduckdb_manifest_error",
                     info = bad)
    }
    for (ok in c("samples", "assay_counts", "a/b/c", "sample_embeddings", "v1..2/x")) {
        expect_silent(.checkResourcePath(ok))
    }
    out <- .safe_package()
    on.exit(unlink(out, recursive = TRUE))
    .edit_manifest(out, function(dp) {
        dp$resources[[1L]]$path <- "../outside/samples"
        dp
    })
    expect_error(readParquet(out), "leaves the package")
})

test_that("a column that the manifest lists and the file lacks is an error, and an unlisted one a warning", {
    out <- .safe_package()
    on.exit(unlink(out, recursive = TRUE))
    sdir <- file.path(out, "samples")
    pq <- list.files(sdir, "[.]parquet$", full.names = TRUE, recursive = TRUE)[[1L]]
    d <- as.data.frame(arrow::read_parquet(pq))
    arrow::write_parquet(cbind(d, extra = 1:2), pq)
    expect_warning(readParquet(out), "does not list, which are not read: 'extra'")
    arrow::write_parquet(d[, setdiff(names(d), "b"), drop = FALSE], pq)
    expect_error(readParquet(out), "has no column for 'b'")
})
