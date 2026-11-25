# Names that become keys, directories or columns are refused when writing would
# rename, merge or misplace them. Each of these used to be written, with the
# name changed or an error from arrow that did not say which name.

skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("SingleCellExperiment")

.entry_se <- function(rownames = c("g1", "g2", "g3"), colnames = c("c1", "c2"),
                      assays = NULL) {
    m <- matrix(c(1, 0, 2, 0, 3, 4), 3, 2)
    SummarizedExperiment::SummarizedExperiment(
        assays = assays %||% list(counts = m),
        rowData = S4Vectors::DataFrame(a = 1:3, row.names = NULL),
        colData = S4Vectors::DataFrame(b = 1:2, row.names = NULL)) |>
        (\(se) { rownames(se) <- rownames; colnames(se) <- colnames; se })()
}

.entry_write <- function(x) {
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE))
    writeParquet(x, out)
}

test_that("repeated, NA and empty dimnames of an experiment are an error", {
    expect_error(.entry_write(.entry_se(rownames = c("g", "g", "h"))),
                 "row names of 'x' have repeats.*'g'")
    expect_error(.entry_write(.entry_se(colnames = c("c", "c"))),
                 "column names of 'x' have repeats")
    expect_error(.entry_write(.entry_se(rownames = c("g", NA, "h"))),
                 "row names of 'x' must not be NA or empty.*positions 2")
    expect_error(.entry_write(.entry_se(colnames = c("c", ""))),
                 "column names of 'x' must not be NA or empty")
})

test_that("an experiment with no dimnames is numbered, and names that differ are kept", {
    se <- .entry_se()
    rownames(se) <- NULL
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE))
    expect_no_error(writeParquet(se, out))
    # g and G are two rows: the names are values, not column names
    out2 <- tempfile()
    on.exit(unlink(out2, recursive = TRUE), add = TRUE)
    expect_no_error(writeParquet(.entry_se(rownames = c("g", "G", "h")), out2))
})

test_that("assay names that cannot be told apart or become paths are an error", {
    m <- matrix(c(1, 0, 2, 0, 3, 4), 3, 2)
    expect_error(.entry_write(.entry_se(assays = list(counts = m, Counts = m))),
                 "assays include names that are the same, ignoring case")
    expect_error(.entry_write(.entry_se(assays = list(`a/b` = m))),
                 "assays become file names")
    expect_error(.entry_write(.entry_se(assays = list(m))), "assays have no names")
    expect_error(.entry_write(.entry_se(assays = list(counts = m, m))),
                 "must not be NA or empty")
})

test_that("reduced dimension names that cannot be told apart are an error", {
    skip_if_not_installed("SingleCellExperiment")
    m <- matrix(c(1, 0, 2, 0, 3, 4), 3, 2, dimnames = list(paste0("g", 1:3), paste0("c", 1:2)))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
    SingleCellExperiment::reducedDims(sce) <- list(PCA = matrix(1:4, 2), pca = matrix(5:8, 2))
    expect_error(.entry_write(sce), "reduced dimensions cannot be written.*'PCA', 'pca'")
})

test_that("a data frame with an empty or repeated column name is an error", {
    df <- S4Vectors::DataFrame(a = 1:2, b = 3:4, check.names = FALSE)
    colnames(df)[2L] <- ""
    expect_error(.entry_write(df), "must not be empty")
    dup <- S4Vectors::DataFrame(a = 1:2, A = 3:4, check.names = FALSE)
    expect_error(.entry_write(dup), "'a', 'A'")
})

test_that("metadata items with names that are the same up to case are an error", {
    se <- .entry_se()
    S4Vectors::metadata(se) <- list(Info = data.frame(x = 1), info = data.frame(y = 2))
    expect_error(.entry_write(se), "metadata of 'x' has items whose names are the same")
    S4Vectors::metadata(se) <- list(`a/b` = data.frame(x = 1))
    expect_error(.entry_write(se), "metadata items become file names")
})
