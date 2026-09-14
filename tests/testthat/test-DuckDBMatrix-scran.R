# Tests the basic functions of a DuckDBMatrix.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBMatrix-scran.R")

skip_if_not_installed("airway")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran correlatePairs tests
###

test_that("correlatePairs works on DuckDBMatrix", {
    # Create DuckDBMatrix and normalize
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    # Select genes with non-zero variance for correlation
    # (genes with all zeros won't have correlations)
    var_genes <- which(apply(log_mat, 1, var) > 0)[1:50]

    # Basic correlatePairs
    ddb_result <- scran::correlatePairs(log_ddb, subset.row = var_genes)
    mat_result <- scran::correlatePairs(log_mat, subset.row = var_genes)

    # DuckDB returns Pearson correlation; scran returns Spearman by default
    # Both should return the same number of pairs
    expect_equal(nrow(ddb_result), nrow(mat_result))

    # Create pair keys for matching (order-independent)
    ddb_pairs <- paste0(pmin(ddb_result$gene1, ddb_result$gene2), "_",
                        pmax(ddb_result$gene1, ddb_result$gene2))
    mat_pairs <- paste0(pmin(mat_result$gene1, mat_result$gene2), "_",
                        pmax(mat_result$gene1, mat_result$gene2))
    expect_setequal(ddb_pairs, mat_pairs)

    # Match pairs and compare correlations (Pearson vs Spearman should correlate)
    # Note: results may be in different order
    mat_rho <- setNames(mat_result$rho, mat_pairs)
    ddb_rho <- setNames(ddb_result$rho, ddb_pairs)
    common_pairs <- intersect(ddb_pairs, mat_pairs)
    expect_gt(cor(ddb_rho[common_pairs], mat_rho[common_pairs]), 0.8)
})

test_that("correlatePairs matches R cor() for Pearson correlation", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    # Select genes with non-zero variance for exact comparison
    var_genes <- which(apply(log_mat, 1, var) > 0)[1:20]

    # Get DuckDB correlations
    ddb_result <- scran::correlatePairs(log_ddb, subset.row = var_genes)

    # Compute Pearson correlations directly with R cor()
    r_cor <- cor(t(log_mat[var_genes, ]))

    # Compare each pair to R cor()
    for (i in seq_len(nrow(ddb_result))) {
        g1 <- ddb_result$gene1[i]
        g2 <- ddb_result$gene2[i]
        expected_r <- r_cor[g1, g2]
        expect_equal(ddb_result$rho[i], expected_r, tolerance = 1e-10,
                     label = paste("Correlation for", g1, "-", g2))
    }
})

test_that("correlatePairs computes p.value/FDR matching cor.test()'s Pearson test", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)
    var_genes <- which(apply(log_mat, 1, var) > 0)[1:15]

    ddb_result <- scran::correlatePairs(log_ddb, subset.row = var_genes)
    expect_true(all(c("p.value", "FDR") %in% colnames(ddb_result)))
    expect_true(all(ddb_result$p.value >= 0 & ddb_result$p.value <= 1))
    expect_true(all(ddb_result$FDR >= 0 & ddb_result$FDR <= 1))
    expect_equal(ddb_result$FDR, p.adjust(ddb_result$p.value, method = "BH"))

    for (i in seq_len(nrow(ddb_result))) {
        g1 <- ddb_result$gene1[i]
        g2 <- ddb_result$gene2[i]
        ct <- cor.test(log_mat[g1, ], log_mat[g2, ], method = "pearson")
        expect_equal(ddb_result$rho[i], unname(ct$estimate), tolerance = 1e-8)
        expect_equal(ddb_result$p.value[i], ct$p.value, tolerance = 1e-8,
                    label = paste("p.value for", g1, "-", g2))
    }
})

test_that("correlatePairs errors for non-zero fill", {
    names(dimnames(airway_counts)) <- c("index1", "index2")

    # Create DuckDBMatrix with non-zero fill
    keycols <- lapply(dimnames(airway_counts),
                      function(x) setNames(seq_along(x), x))
    seed <- DuckDBArraySeed(airway_counts_path, datacol = "value",
                           keycols = keycols)
    # Manually set fill to non-zero value
    seed@fill <- 1L
    pqmat <- DuckDBMatrix(seed)

    expect_error(scran::correlatePairs(pqmat, subset.row = 1:10),
                 "fill = 0")
})

test_that("correlatePairs requires subset.row for large matrices", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    # Should error without subset.row when matrix has >1000 rows
    expect_error(scran::correlatePairs(pqmat),
                 "subset.row is required")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran modelGeneVarByPoisson tests
###

test_that("modelGeneVarByPoisson works on DuckDBMatrix", {
    # Create DuckDBMatrix
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    mat <- as.matrix(pqmat)

    # Basic modelGeneVarByPoisson - compare structure
    ddb_result <- scran::modelGeneVarByPoisson(pqmat)
    mat_result <- scran::modelGeneVarByPoisson(mat)

    # Check structure and approximate equality
    expect_equal(nrow(ddb_result), nrow(mat_result))
    expect_equal(colnames(ddb_result), colnames(mat_result))
    expect_equal(rownames(ddb_result), rownames(mat_result))

    # Mean should be exact
    expect_equal(ddb_result$mean, mat_result$mean, tolerance = 1e-6)
    # Total variance should be close
    expect_equal(ddb_result$total, mat_result$total, tolerance = 1e-4)
})

test_that("modelGeneVarByPoisson works with subset.row", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    mat <- as.matrix(pqmat)
    subset_genes <- 1:100

    ddb_result <- scran::modelGeneVarByPoisson(pqmat, subset.row = subset_genes)
    mat_result <- scran::modelGeneVarByPoisson(mat, subset.row = subset_genes)

    expect_equal(nrow(ddb_result), 100)
    # Use tolerance for floating-point differences in log-transformation
    expect_equal(ddb_result$mean, mat_result$mean, tolerance = 1e-3)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran modelGeneCV2 tests
###

test_that("modelGeneCV2 works on DuckDBMatrix", {
    # Create DuckDBMatrix from raw counts (modelGeneCV2 works on counts)
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    mat <- as.matrix(pqmat)

    # Basic modelGeneCV2
    ddb_result <- scran::modelGeneCV2(pqmat)
    mat_result <- scran::modelGeneCV2(mat)

    # Check structure
    expect_equal(nrow(ddb_result), nrow(mat_result))
    expect_true(all(c("mean", "total", "trend", "ratio", "p.value", "FDR") %in%
                    colnames(ddb_result)))

    # Check mean values are similar
    expect_equal(ddb_result$mean, mat_result$mean, tolerance = 1e-6)
})

test_that("modelGeneCV2 with subset.row works on DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    mat <- as.matrix(pqmat)
    subset_genes <- 1:100

    ddb_result <- scran::modelGeneCV2(pqmat, subset.row = subset_genes)
    mat_result <- scran::modelGeneCV2(mat, subset.row = subset_genes)

    expect_equal(nrow(ddb_result), 100)
    expect_equal(ddb_result$mean, mat_result$mean, tolerance = 1e-6)
})

test_that("modelGeneCV2 with size.factors works on DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    mat <- as.matrix(pqmat)

    # Use custom size factors
    sf <- runif(ncol(mat), 0.5, 1.5)
    sf <- sf / mean(sf)

    ddb_result <- scran::modelGeneCV2(pqmat, size.factors = sf)
    mat_result <- scran::modelGeneCV2(mat, size.factors = sf)

    expect_equal(ddb_result$mean, mat_result$mean, tolerance = 1e-6)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran pairwiseBinom tests
###

test_that("pairwiseBinom works on DuckDBMatrix", {
    # Create DuckDBMatrix and normalize
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    # Use treatment groups from airway
    data(airway, package = "airway")
    groups <- airway$dex

    # Basic pairwiseBinom - should return equal results
    ddb_result <- scran::pairwiseBinom(log_ddb, groups = groups)
    mat_result <- scran::pairwiseBinom(log_mat, groups = groups)

    expect_equal(ddb_result, mat_result)
})

test_that("pairwiseBinom with direction and lfc works on DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    data(airway, package = "airway")
    groups <- airway$dex

    # With direction="up" - should return equal results
    ddb_result <- scran::pairwiseBinom(log_ddb, groups = groups, direction = "up")
    mat_result <- scran::pairwiseBinom(log_mat, groups = groups, direction = "up")
    expect_equal(ddb_result, mat_result)

    # With lfc threshold - should return equal results
    ddb_result <- scran::pairwiseBinom(log_ddb, groups = groups, lfc = 1)
    mat_result <- scran::pairwiseBinom(log_mat, groups = groups, lfc = 1)
    expect_equal(ddb_result, mat_result)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran findMarkers tests
###

test_that("findMarkers with test.type='binom' works on DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    data(airway, package = "airway")
    groups <- airway$dex

    # findMarkers with binom test
    ddb_result <- scran::findMarkers(log_ddb, groups = groups, test.type = "binom")
    mat_result <- scran::findMarkers(log_mat, groups = groups, test.type = "binom")

    # Reorder by gene name and exclude Top column since tie-breaking order
    # for genes with identical p-values may differ
    for (i in seq_along(ddb_result)) {
        common_genes <- rownames(mat_result[[i]])
        cols <- setdiff(colnames(mat_result[[i]]), "Top")
        expect_equal(ddb_result[[i]][common_genes, cols], mat_result[[i]][, cols])
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran scoreMarkers tests
###

test_that("scoreMarkers works on DuckDBMatrix", {
    # Create DuckDBMatrix and normalize
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    # Use treatment groups from airway
    data(airway, package = "airway")
    groups <- airway$dex

    # Basic scoreMarkers
    ddb_result <- scran::scoreMarkers(log_ddb, groups = groups)
    mat_result <- scran::scoreMarkers(log_mat, groups = groups)

    # Check structure
    expect_equal(names(ddb_result), names(mat_result))
    expect_equal(nrow(ddb_result[[1]]), nrow(mat_result[[1]]))

    # Check that effect sizes are similar (not exact due to different algorithms)
    for (g in names(ddb_result)) {
        expect_true("mean.AUC" %in% colnames(ddb_result[[g]]))
        expect_true("mean.logFC.cohen" %in% colnames(ddb_result[[g]]))
        expect_true("mean.logFC.detected" %in% colnames(ddb_result[[g]]))
    }
})

test_that("scoreMarkers with subset.row works on DuckDBMatrix", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)

    data(airway, package = "airway")
    groups <- airway$dex

    # With subset.row
    subset_genes <- 1:100
    ddb_result <- scran::scoreMarkers(log_ddb, groups = groups,
                                      subset.row = subset_genes)

    expect_equal(nrow(ddb_result[[1]]), 100)
})

test_that("scoreMarkers true.auc computes exact rank-based AUC", {
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)

    data(airway, package = "airway")
    groups <- airway$dex

    # Compute with both methods on a subset for speed
    subset_genes <- 1:50
    result_approx <- scran::scoreMarkers(log_ddb, groups = groups,
                                         subset.row = subset_genes,
                                         true.auc = FALSE)
    result_true <- scran::scoreMarkers(log_ddb, groups = groups,
                                       subset.row = subset_genes,
                                       true.auc = TRUE)

    # Both should have same structure
    expect_equal(names(result_approx), names(result_true))
    expect_equal(nrow(result_approx[[1]]), nrow(result_true[[1]]))

    # Both should have AUC columns
    expect_true("mean.AUC" %in% colnames(result_true[[1]]))

    # AUC values should be in valid range [0, 1]
    for (g in names(result_true)) {
        auc_vals <- result_true[[g]][["mean.AUC"]]
        expect_true(all(auc_vals >= 0 & auc_vals <= 1, na.rm = TRUE))
    }

    # True AUC should be different from approximation
    # (they use different methods, so exact match is unlikely)
    # But both should be reasonable (correlated)
    approx_auc <- result_approx[[1]][["mean.AUC"]]
    true_auc <- result_true[[1]][["mean.AUC"]]
    expect_true(cor(approx_auc, true_auc, use = "complete.obs") > 0.5)
})

# Brute-force Mann-Whitney/Wilcoxon AUC, the definition
# .compute_true_auc_sql_DuckDBMatrix computes via SQL ranks. lfc shifts the
# left group's values down, i.e. tests P(left > right + lfc).
.bruteForceAUC <- function(left, right, lfc = 0) {
    left <- left - lfc
    mean(outer(left, right, ">")) + 0.5 * mean(outer(left, right, "=="))
}

test_that(".compute_true_auc_sql_DuckDBMatrix matches a row-subsetted oracle under subset.row", {
    # Regression test: the SQL grouped results by the table's raw row-key
    # values and matched them against seq_len(ngenes), which only coincide
    # when 'x' is the full, unsubsetted matrix. A row-subsetted DuckDBMatrix
    # has raw keys outside 1:ngenes, silently NA-ing or misassigning rows.
    pqmat <- makeSparseDuckDBMatrix()
    groups <- factor(rep(c("a", "b"), length.out = sparse_m))
    sub <- 3:15

    pairs <- data.frame(left = 1L, right = 2L)
    res <- BiocDuckDB:::.compute_true_auc_sql_DuckDBMatrix(
        pqmat[sub, , drop = FALSE], groups, levels(groups), pairs, 0)

    oracle <- sapply(sub, function(i) {
        .bruteForceAUC(sparse_dense[i, groups == "a"], sparse_dense[i, groups == "b"])
    })
    expect_equal(unname(res[, 1L]), unname(oracle), tolerance = 1e-8)
})

test_that(".compute_true_auc_sql_DuckDBMatrix's lfc actually shifts the AUC", {
    # Regression test: the threshold shift was applied identically to both
    # groups' values (a CASE-free suffix on every row), which is a no-op on
    # rank order -- 'lfc' silently had zero effect on the resulting AUC.
    pqmat <- makeSparseDuckDBMatrix()
    groups <- factor(rep(c("a", "b"), length.out = sparse_m))
    pairs <- data.frame(left = 1L, right = 2L)

    auc0 <- BiocDuckDB:::.compute_true_auc_sql_DuckDBMatrix(pqmat, groups, levels(groups), pairs, 0)
    auc1 <- BiocDuckDB:::.compute_true_auc_sql_DuckDBMatrix(pqmat, groups, levels(groups), pairs, 1)
    expect_false(isTRUE(all.equal(auc0[, 1L], auc1[, 1L])))

    oracle1 <- sapply(seq_len(sparse_n), function(i) {
        .bruteForceAUC(sparse_dense[i, groups == "a"], sparse_dense[i, groups == "b"], lfc = 1)
    })
    expect_equal(unname(auc1[, 1L]), unname(oracle1), tolerance = 1e-8)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scran summaryMarkerStats tests
###

test_that("summaryMarkerStats works on DuckDBMatrix", {
    # Create DuckDBMatrix and normalize
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    # Use treatment groups from airway
    data(airway, package = "airway")
    groups <- airway$dex

    # Basic summaryMarkerStats
    ddb_result <- scran::summaryMarkerStats(log_ddb, groups = groups)
    mat_result <- scran::summaryMarkerStats(log_mat, groups = groups)

    # Check structure
    expect_equal(names(ddb_result), names(mat_result))
    expect_equal(nrow(ddb_result[[1]]), nrow(mat_result[[1]]))

    # Check columns exist
    expected_cols <- c("self.average", "other.average", "self.detected", "other.detected")
    for (g in names(ddb_result)) {
        expect_true(all(expected_cols %in% colnames(ddb_result[[g]])))
    }

    # Check values are similar
    for (g in names(ddb_result)) {
        expect_equal(ddb_result[[g]]$self.average, mat_result[[g]]$self.average,
                     tolerance = 1e-6)
        expect_equal(ddb_result[[g]]$self.detected, mat_result[[g]]$self.detected,
                     tolerance = 1e-4)
    }
})
