# Confirms DuckDBIrlbaParam's fast path (materialize once, drive irlba
# directly via BiocSingular:::runIrlbaSVD) matches the same oracle as the
# ordinary lazy BSPARAM = IrlbaParam() path (test-DuckDBMatrix-svd.R), while
# being meaningfully faster, and that it never breaks previously-working
# (if slower) code when it doesn't apply.
#
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBIrlbaParam.R")

skip_if_not_installed("BiocSingular")
skip_if_not_installed("irlba")

test_that("runSVD via DuckDBIrlbaParam matches the dense irlba oracle", {
    pqmat <- makeSparseDuckDBMatrix()
    k <- 3L

    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = k, center = TRUE, scale. = FALSE)
    set.seed(1L)
    result <- BiocSingular::runSVD(pqmat, k = k, center = TRUE, scale = FALSE,
                                   BSPARAM = DuckDBIrlbaParam())

    expect_equal(result$d, oracle$d, tolerance = 1e-8)
})

test_that("scran::fixedPCA via DuckDBIrlbaParam matches the oracle", {
    skip_if_not_installed("scran")
    pqmat <- makeSparseDuckDBMatrix()
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = pqmat))

    set.seed(1L)
    oracle <- irlba::irlba(t(sparse_dense), nv = 3L, center = TRUE)

    set.seed(1L)
    res <- withCallingHandlers(
        scran::fixedPCA(sce, rank = 3L, BSPARAM = DuckDBIrlbaParam()),
        warning = function(w) invokeRestart("muffleWarning")
    )
    pc <- SingleCellExperiment::reducedDim(res, "PCA")
    expect_equal(unname(sqrt(colSums(pc^2))), oracle$d, tolerance = 1e-6)
})

test_that("scater::calculatePCA via DuckDBIrlbaParam matches the oracle", {
    skip_if_not_installed("scater")
    pqmat <- makeSparseDuckDBMatrix()

    set.seed(1L)
    oracle <- irlba::irlba(t(sparse_dense), nv = 3L, center = TRUE)

    set.seed(1L)
    pc <- withCallingHandlers(
        scater::calculatePCA(pqmat, ncomponents = 3L, BSPARAM = DuckDBIrlbaParam()),
        warning = function(w) invokeRestart("muffleWarning")
    )
    expect_equal(unname(sqrt(colSums(pc^2))), oracle$d, tolerance = 1e-6)
})

test_that("DuckDBIrlbaParam falls back gracefully for a non-DuckDBMatrix input", {
    k <- 2L
    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = k, center = TRUE)
    set.seed(1L)
    result <- BiocSingular::runSVD(sparse_dense, k = k, center = TRUE,
                                   BSPARAM = DuckDBIrlbaParam())
    expect_equal(result$d, oracle$d, tolerance = 1e-8)
})

test_that("DuckDBIrlbaParam does not regress a non-zero-filled seed", {
    # A non-zero-filled DuckDBMatrix is not supported by %*% either (see
    # .matmult_DuckDBMatrix_vector's own "must be a zero-filled array" check
    # in DuckDBArray), so the ordinary lazy BSPARAM = IrlbaParam() path
    # already fails on this input; DuckDBIrlbaParam's fallback correctly
    # reproduces that same pre-existing failure rather than papering over it
    # with an incorrect result.
    pqmat <- makeSparseDuckDBMatrix()
    pqmat@seed@fill <- 1
    k <- 2L
    lazy_err <- tryCatch(
        BiocSingular::runSVD(pqmat, k = k, center = TRUE, BSPARAM = BiocSingular::IrlbaParam()),
        error = function(e) conditionMessage(e))
    fast_err <- tryCatch(
        BiocSingular::runSVD(pqmat, k = k, center = TRUE, BSPARAM = DuckDBIrlbaParam()),
        error = function(e) conditionMessage(e))
    expect_identical(fast_err, lazy_err)
})

test_that("DuckDBIrlbaParam falls back gracefully when memory_limit is exceeded", {
    pqmat <- makeSparseDuckDBMatrix()
    k <- 2L
    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = k, center = TRUE)
    set.seed(1L)
    result <- expect_no_error(
        BiocSingular::runSVD(pqmat, k = k, center = TRUE,
                             BSPARAM = DuckDBIrlbaParam(memory_limit = 1))
    )
    expect_equal(result$d, oracle$d, tolerance = 1e-8)
})

test_that("DuckDBIrlbaParam is meaningfully faster than the lazy pushdown path", {
    pqmat <- makeSparseDuckDBMatrix()
    k <- 2L
    set.seed(1L)
    t_fast <- system.time(
        BiocSingular::runSVD(pqmat, k = k, center = TRUE, BSPARAM = DuckDBIrlbaParam())
    )[["elapsed"]]
    set.seed(1L)
    t_lazy <- system.time(
        BiocSingular::runSVD(pqmat, k = k, center = TRUE, BSPARAM = BiocSingular::IrlbaParam())
    )[["elapsed"]]
    expect_lt(t_fast, t_lazy)
})
