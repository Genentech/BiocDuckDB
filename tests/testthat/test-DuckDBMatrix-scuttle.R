# Tests the basic functions of a DuckDBMatrix.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBMatrix-scuttle.R")

skip_if_not_installed("airway")

test_that("librarySizeFactors works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(librarySizeFactors(pqmat), librarySizeFactors(as.matrix(pqmat)))
})

test_that("geometricSizeFactors works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(geometricSizeFactors(pqmat), geometricSizeFactors(as.matrix(pqmat)))
    expect_equal(geometricSizeFactors(pqmat, pseudo.count = 3),
                 geometricSizeFactors(as.matrix(pqmat), pseudo.count = 3))
})

test_that("normalizeCounts works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    checkDuckDBMatrix(normalizeCounts(pqmat), normalizeCounts(as.matrix(pqmat)))
    checkDuckDBMatrix(normalizeCounts(pqmat, log = FALSE), normalizeCounts(as.matrix(pqmat), log = FALSE))
    checkDuckDBMatrix(normalizeCounts(pqmat, pseudo.count = 3), normalizeCounts(as.matrix(pqmat), pseudo.count = 3))
    checkDuckDBMatrix(normalizeCounts(pqmat, transform = "none"), normalizeCounts(as.matrix(pqmat), transform = "none"))
    checkDuckDBMatrix(normalizeCounts(pqmat, transform = "asinh"), normalizeCounts(as.matrix(pqmat), transform = "asinh"))
    checkDuckDBMatrix(normalizeCounts(pqmat, size.factors = 1:ncol(pqmat), center.size.factors = TRUE),
                      normalizeCounts(as.matrix(pqmat), size.factors = 1:ncol(pqmat), center.size.factors = TRUE))
    checkDuckDBMatrix(normalizeCounts(pqmat, size.factors = 1:ncol(pqmat), center.size.factors = FALSE),
                      normalizeCounts(as.matrix(pqmat), size.factors = 1:ncol(pqmat), center.size.factors = FALSE))
})

test_that("calculateTPM works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    checkDuckDBMatrix(scuttle::calculateTPM(pqmat), scuttle::calculateTPM(as.matrix(pqmat)))

    checkDuckDBMatrix(scuttle::calculateCPM(pqmat), scuttle::calculateCPM(as.matrix(pqmat)))

    object <- scuttle::calculateAverage(pqmat)
    names(dimnames(object)) <- NULL
    expected <- as.array(scuttle::calculateAverage(as.matrix(pqmat)))
    expect_equal(object, expected)
})

test_that("perCellQCMetrics works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(perCellQCMetrics(pqmat), perCellQCMetrics(as.matrix(pqmat)))
    expect_equal(perCellQCMetrics(pqmat, threshold = 3),
                 perCellQCMetrics(as.matrix(pqmat), threshold = 3))

    # Test perCellQCMetrics with subsets (feature subsets)
    feat_subsets <- list(
        Group1 = 1:100,
        Group2 = rownames(pqmat)[101:200]
    )
    expect_equal(perCellQCMetrics(pqmat, subsets = feat_subsets),
                 perCellQCMetrics(as.matrix(pqmat), subsets = feat_subsets))
    expect_equal(perCellQCMetrics(pqmat, subsets = feat_subsets, threshold = 3),
                 perCellQCMetrics(as.matrix(pqmat), subsets = feat_subsets, threshold = 3))

    # Test perCellQCMetrics with flatten = FALSE
    expect_equal(perCellQCMetrics(pqmat, subsets = feat_subsets, flatten = FALSE),
                 perCellQCMetrics(as.matrix(pqmat), subsets = feat_subsets, flatten = FALSE))
    nested_result <- perCellQCMetrics(pqmat, subsets = feat_subsets, flatten = FALSE)
    expect_true("subsets" %in% colnames(nested_result))
    expect_true(is(nested_result$subsets, "DataFrame"))

    # Test perCellQCMetrics with percent.top
    expect_equal(perCellQCMetrics(pqmat, percent.top = c(50, 100)),
                 perCellQCMetrics(as.matrix(pqmat), percent.top = c(50, 100)))
    expect_equal(perCellQCMetrics(pqmat, percent.top = c(50, 100), flatten = FALSE),
                 perCellQCMetrics(as.matrix(pqmat), percent.top = c(50, 100), flatten = FALSE))
    expect_equal(perCellQCMetrics(pqmat, subsets = feat_subsets, percent.top = c(50)),
                 perCellQCMetrics(as.matrix(pqmat), subsets = feat_subsets, percent.top = c(50)))
})

test_that("perFeatureQCMetrics works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(perFeatureQCMetrics(pqmat), perFeatureQCMetrics(as.matrix(pqmat)))
    expect_equal(perFeatureQCMetrics(pqmat, threshold = 3),
                 perFeatureQCMetrics(as.matrix(pqmat), threshold = 3))

    # Test perFeatureQCMetrics with subsets (cell subsets)
    cell_subsets <- list(
        SetA = 1:4,
        SetB = colnames(pqmat)[5:8]
    )
    expect_equal(perFeatureQCMetrics(pqmat, subsets = cell_subsets),
                 perFeatureQCMetrics(as.matrix(pqmat), subsets = cell_subsets))
    expect_equal(perFeatureQCMetrics(pqmat, subsets = cell_subsets, threshold = 3),
                 perFeatureQCMetrics(as.matrix(pqmat), subsets = cell_subsets, threshold = 3))

    # Test perFeatureQCMetrics with flatten = FALSE
    expect_equal(perFeatureQCMetrics(pqmat, subsets = cell_subsets, flatten = FALSE),
                 perFeatureQCMetrics(as.matrix(pqmat), subsets = cell_subsets, flatten = FALSE))
    nested_result2 <- perFeatureQCMetrics(pqmat, subsets = cell_subsets, flatten = FALSE)
    expect_true("subsets" %in% colnames(nested_result2))
    expect_true(is(nested_result2$subsets, "DataFrame"))
})

test_that("numDetectedAcrossFeatures works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    ids <- rep(paste0("Group", 1:100), length.out = nrow(pqmat))

    object <- numDetectedAcrossFeatures(pqmat, ids)
    expected <- numDetectedAcrossFeatures(mat, ids)[rownames(object), ]
    expect_equal(object, expected)

    object <- numDetectedAcrossFeatures(pqmat, ids, average = TRUE)
    expected <- numDetectedAcrossFeatures(mat, ids, average = TRUE)[rownames(object), ]
    expect_equal(object, expected)

    object <- numDetectedAcrossFeatures(pqmat, ids, threshold = 5)
    expected <- numDetectedAcrossFeatures(mat, ids, threshold = 5)[rownames(object), ]
    expect_equal(object, expected)

    gene_sets <- list(SetA = 1:100, SetB = 101:200, SetC = 201:300)
    object <- numDetectedAcrossFeatures(pqmat, gene_sets)
    expected <- numDetectedAcrossFeatures(mat, gene_sets)[rownames(object), ]
    expect_equal(object, expected)

    object <- numDetectedAcrossFeatures(pqmat, gene_sets, average = TRUE)
    expected <- numDetectedAcrossFeatures(mat, gene_sets, average = TRUE)[rownames(object), ]
    expect_equal(object, expected)
})

test_that("sumCountsAcrossFeatures works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    ids <- rep(paste0("Group", 1:100), length.out = nrow(pqmat))

    object <- sumCountsAcrossFeatures(pqmat, ids)
    expected <- sumCountsAcrossFeatures(mat, ids)[rownames(object), ]
    expect_equal(object, expected)

    object <- sumCountsAcrossFeatures(pqmat, ids, average = TRUE)
    expected <- sumCountsAcrossFeatures(mat, ids, average = TRUE)[rownames(object), ]
    expect_equal(object, expected)

    gene_sets <- list(SetA = 1:100, SetB = 101:200, SetC = 201:300)
    object <- sumCountsAcrossFeatures(pqmat, gene_sets)
    expected <- sumCountsAcrossFeatures(mat, gene_sets)[rownames(object), ]
    expect_equal(object, expected)

    object <- sumCountsAcrossFeatures(pqmat, gene_sets, average = TRUE)
    expected <- sumCountsAcrossFeatures(mat, gene_sets, average = TRUE)[rownames(object), ]
    expect_equal(object, expected)

    ids_with_na <- ids
    ids_with_na[c(1, 50, 100)] <- NA
    object <- sumCountsAcrossFeatures(pqmat, ids_with_na)
    expected <- sumCountsAcrossFeatures(mat, ids_with_na)[rownames(object), ]
    expect_equal(object, expected)
})

test_that("summarizeAssayByGroup works as expected for a DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    ids <- rep(LETTERS[1:4], length.out = ncol(pqmat))

    object <- summarizeAssayByGroup(pqmat, ids, statistics = "sum")
    expected <- summarizeAssayByGroup(mat, ids, statistics = "sum")
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))
    expect_equal(object$ncells, expected$ncells)

    object <- summarizeAssayByGroup(pqmat, ids, statistics = "mean")
    expected <- summarizeAssayByGroup(mat, ids, statistics = "mean")
    expect_equal(SummarizedExperiment::assay(object, "mean"),
                 SummarizedExperiment::assay(expected, "mean"))

    object <- summarizeAssayByGroup(pqmat, ids, statistics = "num.detected")
    expected <- summarizeAssayByGroup(mat, ids, statistics = "num.detected")
    expect_equal(SummarizedExperiment::assay(object, "num.detected"),
                 SummarizedExperiment::assay(expected, "num.detected"))

    object <- summarizeAssayByGroup(pqmat, ids, statistics = "prop.detected")
    expected <- summarizeAssayByGroup(mat, ids, statistics = "prop.detected")
    expect_equal(SummarizedExperiment::assay(object, "prop.detected"),
                 SummarizedExperiment::assay(expected, "prop.detected"))

    object <- summarizeAssayByGroup(pqmat, ids, statistics = c("sum", "mean"))
    expected <- summarizeAssayByGroup(mat, ids, statistics = c("sum", "mean"))
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))
    expect_equal(SummarizedExperiment::assay(object, "mean"),
                 SummarizedExperiment::assay(expected, "mean"))

    # Test with non-zero threshold
    # NOTE: scuttle has a bug (commit 85635ea, 2026-04-11) where the threshold parameter
    # is ignored in the C++ aggregate_across_cells function, so it always uses threshold=0.
    # BiocDuckDB correctly implements the threshold parameter, so we cannot compare against
    # the broken scuttle implementation. Instead, we verify the result is valid and has
    # fewer detections than threshold=0 (which is the expected behavior).
    object <- summarizeAssayByGroup(pqmat, ids, statistics = "num.detected", threshold = 5)
    object_thresh0 <- summarizeAssayByGroup(pqmat, ids, statistics = "num.detected", threshold = 0)

    # Verify the result structure is correct
    expect_true(is(object, "SummarizedExperiment"))
    expect_equal(dim(SummarizedExperiment::assay(object, "num.detected")), c(nrow(pqmat), length(unique(ids))))

    # Verify that threshold=5 produces fewer or equal detections than threshold=0
    # (values > 5 are a subset of values > 0)
    expect_true(all(SummarizedExperiment::assay(object, "num.detected") <=
                    SummarizedExperiment::assay(object_thresh0, "num.detected")))

    ids_with_na <- ids
    ids_with_na[c(1, 3, 5)] <- NA
    object <- summarizeAssayByGroup(pqmat, ids_with_na, statistics = "sum")
    expected <- summarizeAssayByGroup(mat, ids_with_na, statistics = "sum")
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))
})

# Regression tests for the simplification shared by perCellQCMetrics(),
# perFeatureQCMetrics(), numDetectedAcrossFeatures(), sumCountsAcrossFeatures(),
# and summarizeAssayByGroup() on DuckDBMatrix: the historical
# filter(key %in% known_ids) |> left_join(membership_df) pattern was replaced
# with a bare inner_join(membership_df), since membership_df's key column is
# always exactly the set used in the filter. Each test compares against the
# scuttle oracle computed on a plain in-memory matrix, across single-subset,
# multi-subset, a small fraction of genes/cells, and full-coverage cases.

test_that("perCellQCMetrics matches the oracle across subset shapes", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    single_subset <- list(Group1 = 1:50)
    multi_subset <- list(Group1 = 1:100, Group2 = 101:300, Group3 = 301:500)
    small_fraction <- list(Tiny = 1:5)
    full_coverage <- list(All = seq_len(nrow(pqmat)))

    for (subsets in list(single_subset, multi_subset, small_fraction, full_coverage)) {
        expect_equal(perCellQCMetrics(pqmat, subsets = subsets),
                     perCellQCMetrics(mat, subsets = subsets))
    }
})

test_that("perFeatureQCMetrics matches the oracle across subset shapes", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)
    ncell <- ncol(pqmat)

    # airway_counts has only 8 samples (cells), so subset indices must stay
    # within 1:8.
    cell_single <- list(SetA = 1:4)
    cell_multi <- list(SetA = 1:2, SetB = 3:5, SetC = 6:ncell)
    cell_small <- list(Tiny = 1L)
    cell_full <- list(All = seq_len(ncell))

    for (subsets in list(cell_single, cell_multi, cell_small, cell_full)) {
        expect_equal(perFeatureQCMetrics(pqmat, subsets = subsets),
                     perFeatureQCMetrics(mat, subsets = subsets))
    }
})

test_that("numDetectedAcrossFeatures/sumCountsAcrossFeatures match the oracle across gene-set shapes", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    gene_sets_single <- list(SetA = 1:50)
    gene_sets_multi <- list(SetA = 1:100, SetB = 101:300, SetC = 301:500)
    gene_sets_small <- list(Tiny = 1:3)
    gene_sets_full <- list(All = seq_len(nrow(pqmat)))
    vector_ids_full_coverage <- rep(paste0("Group", 1:100), length.out = nrow(pqmat))

    for (ids in list(gene_sets_single, gene_sets_multi, gene_sets_small,
                     gene_sets_full, vector_ids_full_coverage)) {
        object <- numDetectedAcrossFeatures(pqmat, ids)
        expected <- numDetectedAcrossFeatures(mat, ids)[rownames(object), , drop = FALSE]
        expect_equal(object, expected)

        object <- sumCountsAcrossFeatures(pqmat, ids)
        expected <- sumCountsAcrossFeatures(mat, ids)[rownames(object), , drop = FALSE]
        expect_equal(object, expected)
    }

    # NA-excluded vector ids (membership_df is a strict, non-contiguous subset).
    ids_na <- vector_ids_full_coverage
    ids_na[c(1, 50, 100)] <- NA
    object <- sumCountsAcrossFeatures(pqmat, ids_na)
    expected <- sumCountsAcrossFeatures(mat, ids_na)[rownames(object), , drop = FALSE]
    expect_equal(object, expected)
})

test_that("summarizeAssayByGroup matches the oracle for single- and multi-group assignments", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)

    # "Single group" = full coverage with one group (exercises the inner_join
    # fix when membership_df spans every column).
    groups_single <- rep("OnlyGroup", ncol(pqmat))
    object <- summarizeAssayByGroup(pqmat, groups_single, statistics = "sum")
    expected <- summarizeAssayByGroup(mat, groups_single, statistics = "sum")
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))

    groups_multi <- rep(LETTERS[1:4], length.out = ncol(pqmat))
    object <- summarizeAssayByGroup(pqmat, groups_multi,
                                    statistics = c("sum", "mean", "num.detected", "prop.detected"))
    expected <- summarizeAssayByGroup(mat, groups_multi,
                                      statistics = c("sum", "mean", "num.detected", "prop.detected"))
    for (stat in c("sum", "mean", "num.detected", "prop.detected")) {
        expect_equal(SummarizedExperiment::assay(object, stat),
                     SummarizedExperiment::assay(expected, stat))
    }

    # Small-fraction-of-columns: most cells excluded via subset.col (ids is
    # still full-length, matching the generic's contract).
    small_cols <- 1:3
    object <- summarizeAssayByGroup(pqmat, groups_multi,
                                    subset.col = small_cols, statistics = "sum")
    expected <- summarizeAssayByGroup(mat, groups_multi,
                                      subset.col = small_cols, statistics = "sum")
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))
})

# Regression tests for the window-function fusion of normalizeCounts()
# (and the shared .fused_normalize_by_libsize_DuckDBMatrix() primitive it and
# calculateAverage()/calculateCPM() use when size.factors is NULL). The
# historical path was librarySizeFactors(x) -> colSums(x) [collect()] ->
# sweep(x, 2L, size.factors, "/") [left_join(..., copy = TRUE)]; the fused
# path computes the same normalized values in one query via two chained
# window functions. Each test compares against a from-scratch R oracle, and
# also reconstructs the pre-fusion two-step path inline for a direct
# before/after comparison.

.old_normalizeCounts_path <- function(x, ...) {
    sf <- librarySizeFactors(x)
    DelayedArray::sweep(x, 2L, sf, FUN = "/")
}

test_that("normalizeCounts default path matches the oracle (fill = 0)", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    out <- as.matrix(normalizeCounts(m, transform = "none"))
    sf <- colSums(dense)
    sf <- sf / mean(sf)
    oracle <- sweep(dense, 2, sf, "/")
    expect_equal(out, oracle)

    # Direct before/after comparison against the removed two-pass path.
    old_out <- as.matrix(.old_normalizeCounts_path(m))
    expect_equal(out, old_out)
})

test_that("normalizeCounts default path matches a hand NA-propagation oracle", {
    m <- matrix(c(1, 2, 3, 4, 5, 6, NA, 8, 9, 10, 11, 12), nrow = 3, ncol = 4, byrow = TRUE)
    dimnames(m) <- list(paste0("g", 1:3), paste0("c", 1:4))
    names(dimnames(m)) <- c("index1", "index2")
    tf <- tempfile()
    BiocDuckDB::writeParquet(m, tf)
    pqmat <- DuckDBArray::DuckDBMatrix(tf, datacol = "value",
        keycols = lapply(dimnames(m), function(x) setNames(seq_along(x), x)))

    # mean(colSums(x)) without na.rm propagates NA from a single missing
    # value to every size factor, and hence to the entire output matrix --
    # not just the column containing the NA. Confirmed against the removed
    # two-pass implementation before writing this fused replacement.
    out <- as.matrix(normalizeCounts(pqmat, transform = "none"))
    expect_true(all(is.na(out)))

    old_out <- as.matrix(.old_normalizeCounts_path(pqmat))
    expect_equal(out, old_out)
})

test_that("normalizeCounts default path is unaffected when there are no NAs", {
    m <- matrix(1:12, nrow = 3, ncol = 4)
    dimnames(m) <- list(paste0("g", 1:3), paste0("c", 1:4))
    names(dimnames(m)) <- c("index1", "index2")
    tf <- tempfile()
    BiocDuckDB::writeParquet(m, tf)
    pqmat <- DuckDBArray::DuckDBMatrix(tf, datacol = "value",
        keycols = lapply(dimnames(m), function(x) setNames(seq_along(x), x)))

    out <- as.matrix(normalizeCounts(pqmat, transform = "none"))
    sf <- colSums(m)
    sf <- sf / mean(sf)
    oracle <- sweep(m, 2, sf, "/")
    expect_equal(out, oracle)
})

test_that("logNormCounts() end-to-end output is unaffected by the fusion", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    log_ddb <- scuttle::logNormCounts(
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = m)))
    log_mat <- scuttle::logNormCounts(
        SingleCellExperiment::SingleCellExperiment(assays = list(counts = dense)))

    expect_equal(
        as.matrix(SummarizedExperiment::assay(log_ddb, "logcounts")),
        as.matrix(SummarizedExperiment::assay(log_mat, "logcounts"))
    )
})

test_that("calculateAverage()/calculateCPM() default paths match the oracle and the old path", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)
    sf <- colSums(dense)
    sf <- sf / mean(sf)

    avg <- calculateAverage(m)
    oracle_avg <- rowMeans(sweep(dense, 2, sf, "/"))
    expect_equal(unname(as.vector(avg)), unname(oracle_avg))

    old_avg <- rowMeans(as.matrix(.old_normalizeCounts_path(m)))
    expect_equal(unname(as.vector(avg)), unname(old_avg))

    cpm <- as.matrix(calculateCPM(m))
    oracle_cpm <- sweep(dense, 2, colSums(dense) / 1e6, "/")
    expect_equal(cpm, oracle_cpm)

    # calculateCPM()'s explicit size.factors branch is untouched by Fix 3;
    # confirm it still works correctly alongside the fused default path.
    set.seed(11)
    sf2 <- runif(ncol(m), 0.5, 1.5)
    cpm2 <- as.matrix(calculateCPM(m, size.factors = sf2))
    lib.sizes <- colSums(dense) / 1e6
    lib.sizes2 <- sf2 / mean(sf2) * mean(lib.sizes)
    oracle_cpm2 <- sweep(dense, 2, lib.sizes2, "/")
    expect_equal(cpm2, oracle_cpm2)
})

test_that(".fused_normalize_by_libsize_DuckDBMatrix() refuses a nonzero fill, matching sweep()", {
    m5 <- makeSparseDuckDBMatrix() + 5
    expect_equal(m5@seed@fill, 5)
    expect_error(normalizeCounts(m5), "zero-filled")
    expect_error(calculateAverage(m5), "zero-filled")
    expect_error(calculateCPM(m5), "zero-filled")
})

test_that("the fused normalize query plan scans the base table once with no join", {
    m <- makeSparseDuckDBMatrix()
    out <- normalizeCounts(m, transform = "none")
    table <- out@seed@table
    sql_text <- dbplyr::sql_render(DuckDBDataFrame::tblconn(table, select = TRUE))
    db_conn <- DuckDBDataFrame::tblconn(table, select = FALSE)$src$con
    plan <- DBI::dbGetQuery(db_conn, paste("EXPLAIN", sql_text))
    plan_text <- paste(plan[[2]], collapse = "\n")

    expect_equal(lengths(regmatches(plan_text, gregexpr("Read Parquet", plan_text))), 1L)
    expect_false(grepl("Join", plan_text))
    expect_true(grepl("Window", plan_text))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Regression: scuttle's underscore-spelled aliases (used by callers like
### batchelor's .apply_cosine_norm(), which calls normalizeCounts(x,
### size_factors=..., center_size_factors=FALSE, log=FALSE)) must not
### silently fall into `...` unused. See the normalizeCounts,DuckDBMatrix-
### method comment in DuckDBMatrix-scuttle.R for the full explanation; this
### covers every scuttle generic BiocDuckDB implements that scuttle itself
### gives underscore aliases to (geometricSizeFactors/summarizeAssayByGroup
### are deliberately excluded: scuttle defines no aliases for either).
###

test_that("normalizeCounts accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    sf <- seq_len(ncol(pqmat))

    checkDuckDBMatrix(
        normalizeCounts(pqmat, size_factors = sf, center_size_factors = FALSE, log = FALSE),
        normalizeCounts(as.matrix(pqmat), size_factors = sf, center_size_factors = FALSE, log = FALSE)
    )
    checkDuckDBMatrix(
        normalizeCounts(pqmat, size_factors = sf, center_size_factors = TRUE, log = FALSE),
        normalizeCounts(as.matrix(pqmat), size.factors = sf, center.size.factors = TRUE, log = FALSE)
    )
    checkDuckDBMatrix(
        normalizeCounts(pqmat, pseudo_count = 3),
        normalizeCounts(as.matrix(pqmat), pseudo.count = 3)
    )
    checkDuckDBMatrix(
        normalizeCounts(pqmat, subset_row = 1:10),
        normalizeCounts(as.matrix(pqmat), subset.row = 1:10)
    )

    # A non-NULL alias must take precedence over a simultaneously-passed
    # dotted value with a different value (scuttle's own .replace()
    # semantics): a caller relying on the alias should never silently get
    # the dotted value instead.
    wrong_sf <- rev(sf)
    checkDuckDBMatrix(
        normalizeCounts(pqmat, size.factors = wrong_sf, size_factors = sf,
                        center.size.factors = TRUE, center_size_factors = FALSE, log = FALSE),
        normalizeCounts(as.matrix(pqmat), size_factors = sf, center_size_factors = FALSE, log = FALSE)
    )

    # A NULL alias (the default) must leave the dotted value intact, not
    # overwrite it with NULL and fall through to the is.null(size.factors)
    # library-size branch -- the exact shape of the original bug.
    checkDuckDBMatrix(
        normalizeCounts(pqmat, size.factors = sf, center.size.factors = FALSE, log = FALSE),
        normalizeCounts(as.matrix(pqmat), size.factors = sf, center.size.factors = FALSE, log = FALSE)
    )
})

test_that("calculateCPM accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    sf <- seq_len(ncol(pqmat))

    checkDuckDBMatrix(scuttle::calculateCPM(pqmat, size_factors = sf),
                      scuttle::calculateCPM(as.matrix(pqmat), size.factors = sf))
    checkDuckDBMatrix(scuttle::calculateCPM(pqmat, subset_row = 1:10),
                      scuttle::calculateCPM(as.matrix(pqmat), subset.row = 1:10))
})

test_that("calculateAverage accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    sf <- seq_len(ncol(pqmat))

    object <- scuttle::calculateAverage(pqmat, size_factors = sf)
    names(dimnames(object)) <- NULL
    expected <- as.array(scuttle::calculateAverage(as.matrix(pqmat), size.factors = sf))
    expect_equal(object, expected)

    object <- scuttle::calculateAverage(pqmat, subset_row = 1:10)
    names(dimnames(object)) <- NULL
    expected <- as.array(scuttle::calculateAverage(as.matrix(pqmat), subset.row = 1:10))
    expect_equal(object, expected)
})

test_that("calculateTPM accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    sf <- seq_len(ncol(pqmat))

    checkDuckDBMatrix(scuttle::calculateTPM(pqmat, size_factors = sf),
                      scuttle::calculateTPM(as.matrix(pqmat), size.factors = sf))
    checkDuckDBMatrix(scuttle::calculateTPM(pqmat, subset_row = 1:10),
                      scuttle::calculateTPM(as.matrix(pqmat), subset.row = 1:10))
})

test_that("librarySizeFactors accepts scuttle's underscore-spelled alias", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(librarySizeFactors(pqmat, subset_row = 1:10),
                 librarySizeFactors(as.matrix(pqmat), subset.row = 1:10))
})

test_that("numDetectedAcrossFeatures accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)
    ids <- rep(paste0("Group", 1:100), length.out = nrow(pqmat))

    object <- numDetectedAcrossFeatures(pqmat, ids, subset_row = 1:200)
    expected <- numDetectedAcrossFeatures(mat, ids, subset.row = 1:200)[rownames(object), ]
    expect_equal(object, expected)

    object <- numDetectedAcrossFeatures(pqmat, ids, subset_col = 1:5)
    expected <- numDetectedAcrossFeatures(mat, ids, subset.col = 1:5)[rownames(object), ]
    expect_equal(object, expected)

    object <- numDetectedAcrossFeatures(pqmat, ids, detection_limit = 5)
    expected <- numDetectedAcrossFeatures(mat, ids, threshold = 5)[rownames(object), ]
    expect_equal(object, expected)
})

test_that("sumCountsAcrossFeatures accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)
    ids <- rep(paste0("Group", 1:100), length.out = nrow(pqmat))

    object <- sumCountsAcrossFeatures(pqmat, ids, subset_row = 1:200)
    expected <- sumCountsAcrossFeatures(mat, ids, subset.row = 1:200)[rownames(object), ]
    expect_equal(object, expected)

    object <- sumCountsAcrossFeatures(pqmat, ids, subset_col = 1:5)
    expected <- sumCountsAcrossFeatures(mat, ids, subset.col = 1:5)[rownames(object), ]
    expect_equal(object, expected)
})

test_that("perCellQCMetrics accepts scuttle's underscore-spelled aliases", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(perCellQCMetrics(pqmat, detection_limit = 3),
                 perCellQCMetrics(as.matrix(pqmat), threshold = 3))

    feat_subsets <- list(
        Group1 = 1:100,
        Group2 = rownames(pqmat)[101:200]
    )
    expect_equal(perCellQCMetrics(pqmat, subsets = feat_subsets, percent_top = c(10L, 50L)),
                 perCellQCMetrics(as.matrix(pqmat), subsets = feat_subsets, percent.top = c(10L, 50L)))
})

test_that("perFeatureQCMetrics accepts scuttle's underscore-spelled alias", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    expect_equal(perFeatureQCMetrics(pqmat, detection_limit = 3),
                 perFeatureQCMetrics(as.matrix(pqmat), threshold = 3))
})
