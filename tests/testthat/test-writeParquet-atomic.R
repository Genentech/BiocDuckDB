# A package is written beside its target and moved into place, so that a write
# that fails leaves the target as it was, and an existing package is replaced
# only when it is asked and only when it is a package.

skip_if_not_installed("SummarizedExperiment")

.atomic_se <- function(value = 1, assays = NULL) {
    m <- matrix(c(value, 0, 2, 0, 3, 4), 3, 2,
                dimnames = list(c("g1", "g2", "g3"), c("c1", "c2")))
    SummarizedExperiment::SummarizedExperiment(
        assays = assays %||% list(counts = m),
        colData = S4Vectors::DataFrame(b = 1:2, row.names = c("c1", "c2")))
}

.atomic_counts <- function(path) {
    as.matrix(SummarizedExperiment::assay(readParquet(path), "counts"))
}

test_that("a package is written, and nothing else is left beside it", {
    parent <- tempfile(); dir.create(parent)
    on.exit(unlink(parent, recursive = TRUE))
    out <- file.path(parent, "pkg")
    writeParquet(.atomic_se(), out)
    expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), "pkg")
    expect_true(file.exists(file.path(out, "datapackage.json")))
    expect_equal(unname(.atomic_counts(out)["g1", "c1"]), 1)
    # the parent of the target is made, as it was
    deep <- file.path(parent, "a", "b", "pkg")
    writeParquet(.atomic_se(), deep)
    expect_true(file.exists(file.path(deep, "datapackage.json")))
})

test_that("a package that exists is not replaced unless it is asked", {
    parent <- tempfile(); dir.create(parent)
    on.exit(unlink(parent, recursive = TRUE))
    out <- file.path(parent, "pkg")
    writeParquet(.atomic_se(1), out)
    expect_error(writeParquet(.atomic_se(7), out), "already exists and is not empty.*overwrite = TRUE")
    expect_equal(unname(.atomic_counts(out)["g1", "c1"]), 1)
    writeParquet(.atomic_se(7), out, overwrite = TRUE)
    expect_equal(unname(.atomic_counts(out)["g1", "c1"]), 7)
    expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), "pkg")
    expect_error(writeParquet(.atomic_se(), out, overwrite = NA), "must be TRUE or FALSE")
})

test_that("a write that fails leaves the target as it was, and nothing beside it", {
    parent <- tempfile(); dir.create(parent)
    on.exit(unlink(parent, recursive = TRUE))
    out <- file.path(parent, "pkg")
    writeParquet(.atomic_se(1), out)
    before <- list.files(out, recursive = TRUE)
    m <- matrix(1:6, 3, 2, dimnames = list(c("g1", "g2", "g3"), c("c1", "c2")))
    bad <- .atomic_se(assays = list(counts = m, Counts = m))
    expect_error(writeParquet(bad, out, overwrite = TRUE), "same, ignoring case")
    expect_identical(list.files(out, recursive = TRUE), before)
    expect_equal(unname(.atomic_counts(out)["g1", "c1"]), 1)
    expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), "pkg")
    # a failure on a new path leaves no directory at all
    fresh <- file.path(parent, "fresh")
    expect_error(writeParquet(bad, fresh), "same, ignoring case")
    expect_false(file.exists(fresh))
    expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), "pkg")
    # and the next write is not taken for a nested one
    expect_no_error(writeParquet(.atomic_se(3), fresh))
})

test_that("only a package is replaced, and only when it holds nothing that the manifest does not list", {
    parent <- tempfile(); dir.create(parent)
    on.exit(unlink(parent, recursive = TRUE))
    stray <- file.path(parent, "mydata")
    dir.create(stray)
    writeLines("precious", file.path(stray, "notes.txt"))
    expect_error(writeParquet(.atomic_se(), stray, overwrite = TRUE),
                 "no datapackage.json, so it is not a package")
    expect_identical(readLines(file.path(stray, "notes.txt")), "precious")

    out <- file.path(parent, "pkg")
    writeParquet(.atomic_se(), out)
    writeLines("mine", file.path(out, "my-notes.txt"))
    expect_error(writeParquet(.atomic_se(), out, overwrite = TRUE),
                 "entries that its manifest does not list: 'my-notes.txt'")
    expect_identical(readLines(file.path(out, "my-notes.txt")), "mine")

    # an empty directory is written into without overwrite, and a file is refused
    empty <- file.path(parent, "empty")
    dir.create(empty)
    expect_no_error(writeParquet(.atomic_se(), empty))
    afile <- file.path(parent, "afile")
    writeLines("x", afile)
    expect_error(writeParquet(.atomic_se(), afile, overwrite = TRUE), "not a directory")
})

test_that("the experiments of a MultiAssayExperiment are written in place, inside the staged package", {
    skip_if_not_installed("MultiAssayExperiment")
    parent <- tempfile(); dir.create(parent)
    on.exit(unlink(parent, recursive = TRUE))
    m <- matrix(c(1, 0, 2, 0, 3, 4), 3, 2, dimnames = list(c("g1", "g2", "g3"), c("c1", "c2")))
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(rna = m),
        colData = S4Vectors::DataFrame(row.names = c("c1", "c2")))
    out <- file.path(parent, "mae")
    writeParquet(mae, out)
    expect_identical(list.files(parent, all.files = TRUE, no.. = TRUE), "mae")
    expect_true(file.exists(file.path(out, "datapackage.json")))
    expect_s4_class(readParquet(out), "MultiAssayExperiment")
})
