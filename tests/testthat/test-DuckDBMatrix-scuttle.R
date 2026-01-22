# Tests the basic functions of a DuckDBMatrix.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBMatrix-scuttle.R")

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

    object <- summarizeAssayByGroup(pqmat, ids, statistics = "num.detected", threshold = 5)
    expected <- summarizeAssayByGroup(mat, ids, statistics = "num.detected", threshold = 5)
    expect_equal(SummarizedExperiment::assay(object, "num.detected"),
                 SummarizedExperiment::assay(expected, "num.detected"))

    ids_with_na <- ids
    ids_with_na[c(1, 3, 5)] <- NA
    object <- summarizeAssayByGroup(pqmat, ids_with_na, statistics = "sum")
    expected <- summarizeAssayByGroup(mat, ids_with_na, statistics = "sum")
    expect_equal(SummarizedExperiment::assay(object, "sum"),
                 SummarizedExperiment::assay(expected, "sum"))
})
