# Tests the beachmat/tatami fast path for DuckDBArraySeed (initializeCpp.R,
# initializeOptions.R, loadIntoMemory.R). This is a performance optimization
# for anything that already goes through beachmat::initializeCpp (e.g.
# BiocSingular's compute_center/compute_scale); it must never change results,
# and it must never turn previously-working (if slower) code into an error --
# every unsupported case should fall back to beachmat's generic path.
#
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-initializeCpp.R")

skip_if_not_installed("beachmat")
skip_if_not_installed("BiocSingular")
skip_if_not_installed("irlba")

test_that("initializeCpp on a sparse DuckDBArraySeed returns a real tatami pointer", {
    pqmat <- makeSparseDuckDBMatrix()
    ptr <- beachmat::initializeCpp(pqmat@seed)
    expect_true(inherits(ptr, "externalptr"))
})

test_that("initializeCpp's fast path does not change compute_center/runSVD results", {
    pqmat <- makeSparseDuckDBMatrix()
    # See test-DuckDBMatrix-svd.R: this fixture has near-degenerate singular
    # values, so irlba's unseeded random start vector must be pinned for a
    # reproducible comparison regardless of what ran earlier in the session.
    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = 3L, center = TRUE)
    set.seed(1L)
    res <- BiocSingular::runSVD(pqmat, k = 3L, center = TRUE, scale = FALSE,
                                BSPARAM = BiocSingular::IrlbaParam())
    expect_equal(res$d, oracle$d, tolerance = 1e-8)
})

test_that("loadIntoMemory() works directly on a DuckDBMatrix or its seed", {
    pqmat <- makeSparseDuckDBMatrix()
    expect_true(inherits(BiocDuckDB::loadIntoMemory(pqmat), "externalptr"))
    expect_true(inherits(BiocDuckDB::loadIntoMemory(pqmat@seed), "externalptr"))
})

test_that("a non-zero-filled seed falls back gracefully instead of erroring", {
    pqmat <- makeSparseDuckDBMatrix()
    bad_seed <- pqmat@seed
    bad_seed@fill <- 1
    expect_no_error(ptr <- beachmat::initializeCpp(bad_seed))
    expect_true(inherits(ptr, "externalptr"))
})

test_that("duckdb.realize = FALSE falls back to beachmat's generic path", {
    pqmat <- makeSparseDuckDBMatrix()
    expect_message(
        ptr <- beachmat::initializeCpp(pqmat@seed, duckdb.realize = FALSE),
        "unknown matrix fallback"
    )
    expect_true(inherits(ptr, "externalptr"))
})

test_that("exceeding memory_limit falls back gracefully instead of erroring", {
    pqmat <- makeSparseDuckDBMatrix()
    old <- initializeOptions("memory_limit", 1)  # guaranteed to be exceeded
    on.exit(initializeOptions("memory_limit", old))
    expect_no_error(ptr <- beachmat::initializeCpp(pqmat@seed))
    expect_true(inherits(ptr, "externalptr"))
})

test_that("initializeOptions get/set round-trips", {
    old <- initializeOptions("realize")
    on.exit(initializeOptions("realize", old))
    initializeOptions("realize", FALSE)
    expect_identical(initializeOptions("realize"), FALSE)
})
