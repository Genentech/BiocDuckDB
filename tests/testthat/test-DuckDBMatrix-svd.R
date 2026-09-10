# Confirms BiocDuckDB's existing %*%/crossprod SQL pushdown for DuckDBMatrix
# (DuckDBMatrix-utils.R) already composes correctly with the wider
# Bioconductor PCA ecosystem (BiocSingular::runSVD, scran::fixedPCA,
# scater::runPCA), closing the gap documented in
# docs/reference/algorithm-catalog.md ("BiocDuckDB: storage/accessor ONLY;
# confirmed ZERO runPCA/eigen/svd calls anywhere in the suite"). No new
# production code was needed: BiocSingular::IrlbaParam()'s default
# (fold = Inf) never routes through DuckDBTable's sweep()-based centering (the
# one path that would be numerically wrong for a COO/sparse table, since an
# absent cell's "-mean" contribution would silently vanish); instead it
# computes center/scale via a safe generic reduction and lets irlba's own
# algebraic no-densify centering trick drive %*% directly, which dispatches
# straight to the existing SQL pushdown.
#
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBMatrix-svd.R")

skip_if_not_installed("BiocSingular")
skip_if_not_installed("irlba")

# Genes are rows, samples are columns, and 88% of (gene, sample) cells are
# implicit zero (COO-encoded, not just numerically small) -- exactly the
# shape a naive sweep()-then-aggregate centering would get wrong.
align_sign <- function(mat, ref) {
    flips <- vapply(seq_len(ncol(mat)), function(j) {
        sign(sum(mat[, j] * ref[, j]))
    }, numeric(1))
    flips[flips == 0] <- 1
    sweep(mat, 2, flips, `*`)
}

test_that("BiocSingular::runSVD on a sparse DuckDBMatrix matches the dense irlba oracle", {
    pqmat <- makeSparseDuckDBMatrix()
    k <- 3L

    # irlba draws a random starting vector from R's global RNG when none is
    # supplied; this fixture has near-degenerate singular values (see below),
    # so an unseeded call is sensitive to whatever ran earlier in the test
    # session. Seed each irlba-consuming call independently for determinism.
    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = k, center = TRUE, scale. = FALSE)
    set.seed(1L)
    result <- BiocSingular::runSVD(pqmat, k = k, center = TRUE, scale = FALSE,
                                   BSPARAM = BiocSingular::IrlbaParam())

    expect_equal(result$d, oracle$d, tolerance = 1e-8)
    expect_equal(unname(align_sign(result$u, oracle$u)), unname(oracle$u), tolerance = 1e-6)
    expect_equal(unname(align_sign(result$v, oracle$v)), unname(oracle$v), tolerance = 1e-6)
})

test_that("BiocSingular::runSVD centering is correct even when most cells are implicit zero", {
    # A regression guard for the specific numerical trap: if centering were
    # done by sweeping only the *present* (non-zero) rows of the COO table
    # and then summing, every implicit-zero cell's "-mean" contribution would
    # be silently dropped and the singular values would NOT match a dense
    # centered SVD computed on the true, fully-materialized matrix.
    pqmat <- makeSparseDuckDBMatrix()
    frac_zero <- mean(sparse_dense == 0)
    expect_gt(frac_zero, 0.75)  # fixture is genuinely sparse, not incidentally so

    set.seed(1L)
    oracle <- irlba::irlba(sparse_dense, nv = 2L, center = TRUE, scale. = FALSE)
    set.seed(1L)
    result <- BiocSingular::runSVD(pqmat, k = 2L, center = TRUE, scale = FALSE,
                                   BSPARAM = BiocSingular::IrlbaParam())
    expect_equal(result$d, oracle$d, tolerance = 1e-8)
})

test_that("scran::fixedPCA on a SingleCellExperiment with a DuckDBMatrix assay matches the oracle", {
    skip_if_not_installed("scran")
    pqmat <- makeSparseDuckDBMatrix()
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = pqmat))

    # fixedPCA transposes the assay before centering (centers per-gene, i.e.
    # per COLUMN of t(x)), so the matching oracle centers t(sparse_dense), not
    # sparse_dense itself.
    set.seed(1L)
    oracle <- irlba::irlba(t(sparse_dense), nv = 3L, center = TRUE)

    set.seed(1L)
    res <- withCallingHandlers(
        scran::fixedPCA(sce, rank = 3L, BSPARAM = BiocSingular::IrlbaParam()),
        warning = function(w) invokeRestart("muffleWarning")  # fixedPCA is .Deprecated
    )
    pc <- SingleCellExperiment::reducedDim(res, "PCA")
    expect_equal(dim(pc), c(sparse_m, 3L))
    expect_equal(unname(sqrt(colSums(pc^2))), oracle$d, tolerance = 1e-6)
})

test_that("scater::runPCA on a SingleCellExperiment with a DuckDBMatrix assay matches the oracle", {
    skip_if_not_installed("scater")
    pqmat <- makeSparseDuckDBMatrix()
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = pqmat))
    set.seed(1L)
    oracle <- irlba::irlba(t(sparse_dense), nv = 3L, center = TRUE)

    set.seed(1L)
    res <- withCallingHandlers(
        scater::runPCA(sce, ncomponents = 3L, BSPARAM = BiocSingular::IrlbaParam()),
        warning = function(w) invokeRestart("muffleWarning")
    )
    pc <- SingleCellExperiment::reducedDim(res, "PCA")
    expect_equal(unname(sqrt(colSums(pc^2))), oracle$d, tolerance = 1e-6)
})

test_that("a computed PCA round-trips through writeParquet/readParquet as DuckDBEmbeddings", {
    skip_if_not_installed("scran")
    pqmat <- makeSparseDuckDBMatrix()
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = pqmat))
    set.seed(1L)
    res <- withCallingHandlers(
        scran::fixedPCA(sce, rank = 3L, BSPARAM = BiocSingular::IrlbaParam()),
        warning = function(w) invokeRestart("muffleWarning")
    )
    pc_before <- SingleCellExperiment::reducedDim(res, "PCA")

    outdir <- file.path(tempfile())
    BiocDuckDB::writeParquet(res, outdir)
    back <- BiocDuckDB::readParquet(outdir)
    rd <- SingleCellExperiment::reducedDim(back, "PCA")

    expect_s4_class(rd, "DuckDBEmbeddings")
    m <- as.matrix(rd)[rownames(pc_before), , drop = FALSE]
    expect_equal(unclass(m), unclass(pc_before), tolerance = 1e-8, check.attributes = FALSE)
})
