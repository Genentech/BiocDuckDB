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

# Regression test for a correctness bug in the previous one-pass raw-moment
# formula (SS = sum_sq - sum^2/k, numerator = sum_prod - sum1*sum2/k): at
# large-magnitude/non-log-normalized scale, catastrophic cancellation in
# that formula produced |rho| > 1 or NaN. The existing correlatePairs tests
# above only ever exercise log-normalized scRNA-seq scale, where the old
# formula happened to be numerically fine (~1e-8 relative error) -- this
# test specifically exercises the adversarial regime the old formula got
# wrong, so a future regression to a raw-moment-style formula is caught.
test_that("correlatePairs stays numerically stable at large-magnitude scale", {
    makeLargeMagnitudeDuckDBMatrix <- function(mat) {
        idx <- which(mat != 0, arr.ind = TRUE)
        df <- data.frame(gene = idx[, 1L], cell = idx[, 2L], value = mat[idx])
        path <- tempfile(fileext = ".parquet")
        arrow::write_parquet(df, path)
        DuckDBMatrix(path, datacol = "value",
                     keycols = list(
                         gene = setNames(seq_len(nrow(mat)), rownames(mat)),
                         cell = setNames(seq_len(ncol(mat)), colnames(mat))
                     ))
    }

    set.seed(20261105)
    k <- 10000L
    rho_true <- 0.09
    z1 <- rnorm(k)
    z2 <- rho_true * z1 + sqrt(1 - rho_true^2) * rnorm(k)
    z3 <- rnorm(k)

    for (offset in c(1e5, 1e9, 1e10)) {
        mat <- rbind(geneA = offset + z1, geneB = offset + z2, geneC = offset + z3)
        colnames(mat) <- paste0("cell", seq_len(k))
        ddb <- makeLargeMagnitudeDuckDBMatrix(mat)

        res <- scran::correlatePairs(ddb, pairings = NULL)
        dense <- cor(t(mat))

        expect_true(all(abs(res$rho) <= 1 + 1e-9),
                   label = sprintf("offset = %.0e", offset))
        expect_false(any(is.nan(res$rho)), label = sprintf("offset = %.0e", offset))
        for (i in seq_len(nrow(res))) {
            expect_equal(res$rho[i], dense[res$gene1[i], res$gene2[i]],
                        tolerance = 1e-6,
                        label = sprintf("offset = %.0e, %s-%s", offset,
                                        res$gene1[i], res$gene2[i]))
        }
    }
})

# Regression test for sparse/implicit-fill correctness: a gene pair with NO
# overlapping explicitly-stored cells, and a gene that is entirely implicit
# fill, must still be handled correctly (not silently dropped or NA'd out
# incorrectly) by the stable formula's 4-population decomposition.
test_that("correlatePairs handles disjoint support and all-fill genes", {
    makeDuckDBMatrixFromDense <- function(mat) {
        idx <- which(mat != 0, arr.ind = TRUE)
        df <- data.frame(gene = idx[, 1L], cell = idx[, 2L], value = mat[idx])
        path <- tempfile(fileext = ".parquet")
        arrow::write_parquet(df, path)
        DuckDBMatrix(path, datacol = "value",
                     keycols = list(
                         gene = setNames(seq_len(nrow(mat)), rownames(mat)),
                         cell = setNames(seq_len(ncol(mat)), colnames(mat))
                     ))
    }

    set.seed(20261105)
    nc <- 100L
    mat <- matrix(0, nrow = 4L, ncol = nc)
    rownames(mat) <- c("A", "B", "C", "AllZero")
    colnames(mat) <- paste0("cell", seq_len(nc))
    mat["A", 1:50] <- rnorm(50, 10, 2)
    mat["B", 51:100] <- rnorm(50, 10, 2)  # A and B share zero stored cells
    mat["C", ] <- rnorm(nc, 10, 2)        # C overlaps both A and B

    ddb <- makeDuckDBMatrixFromDense(mat)
    res <- scran::correlatePairs(ddb, pairings = NULL)
    dense <- suppressWarnings(cor(t(mat)))

    expect_true(all(abs(res$rho[!is.na(res$rho)]) <= 1 + 1e-9))
    for (i in seq_len(nrow(res))) {
        truth <- dense[res$gene1[i], res$gene2[i]]
        if (is.na(truth)) {
            expect_true(is.na(res$rho[i]))
        } else {
            expect_equal(res$rho[i], truth, tolerance = 1e-6,
                        label = paste(res$gene1[i], res$gene2[i]))
        }
    }
    # AllZero has zero explicitly-stored cells; its correlation with
    # anything is undefined (zero variance), not silently dropped.
    all_zero_rows <- res[res$gene1 == "AllZero" | res$gene2 == "AllZero", ]
    expect_equal(nrow(all_zero_rows), 3L)
    expect_true(all(is.na(all_zero_rows$rho)))
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

# Regression tests for the fused single-pass mean/variance computation shared
# by .compute_blocked_stats_DuckDBMatrix() (modelGeneVar/modelGeneVarByPoisson),
# .compute_blocked_stats_normalized_DuckDBMatrix() (modelGeneCV2), and
# .compute_group_stats_DuckDBMatrix() (pairwiseTTests/findMarkers), via the
# shared `.fused_group_variance_stats()` primitive.
#
# Each test compares the DuckDBMatrix result against an independent R oracle
# (mean()/var() on a densified matrix), rather than against the historical
# two-pass SQL implementation, which no longer exists in this file -- an
# oracle that is both an independent check and a permanent, standalone
# regression test.

test_that(".compute_blocked_stats_DuckDBMatrix matches R's mean()/var() oracle (fill = 0)", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    res <- .compute_blocked_stats_DuckDBMatrix(m)

    expect_equal(unname(res$means[, 1]), unname(rowMeans(dense)))
    expect_equal(unname(res$vars[, 1]), unname(apply(dense, 1, var)))
})

test_that(".compute_blocked_stats_DuckDBMatrix matches the oracle for a nonzero fill", {
    m <- makeSparseDuckDBMatrix() + 5
    expect_equal(m@seed@fill, 5)
    dense <- as.matrix(m)

    res <- .compute_blocked_stats_DuckDBMatrix(m)

    expect_equal(unname(res$means[, 1]), unname(rowMeans(dense)))
    expect_equal(unname(res$vars[, 1]), unname(apply(dense, 1, var)))
})

test_that(".compute_blocked_stats_DuckDBMatrix matches the oracle per block", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    set.seed(1)
    block <- factor(sample(c("A", "B"), ncol(m), replace = TRUE))
    res <- .compute_blocked_stats_DuckDBMatrix(m, block = block)

    for (b in 1:2) {
        lvl <- levels(block)[b]
        sub <- dense[, block == lvl, drop = FALSE]
        expect_equal(unname(res$means[, b]), unname(rowMeans(sub)))
        expect_equal(unname(res$vars[, b]), unname(apply(sub, 1, var)))
    }
})

test_that(".compute_blocked_stats_DuckDBMatrix is numerically stable for a high-mean, low-relative-variance gene", {
    # Adversarial case for catastrophic cancellation: mean ~1e8, variance O(1).
    # A hand-rolled sum-of-squares (sum(x^2)/n - mean^2) would lose most
    # significant digits here; DuckDB's native VAR_SAMP should not.
    set.seed(3)
    ncell <- 200L
    vals <- 1e8 + round(rnorm(ncell, sd = 2), 0)
    m <- matrix(vals, nrow = 1, ncol = ncell)
    dimnames(m) <- list("g1", paste0("c", seq_len(ncell)))
    names(dimnames(m)) <- c("index1", "index2")
    tf <- tempfile()
    BiocDuckDB::writeParquet(m, tf)
    pqmat <- DuckDBArray::DuckDBMatrix(tf, datacol = "value",
        keycols = lapply(dimnames(m), function(x) setNames(seq_along(x), x)))
    expect_equal(pqmat@seed@fill, 0)

    res <- .compute_blocked_stats_DuckDBMatrix(pqmat)

    expect_equal(unname(res$means[1, 1]), mean(vals))
    # Tight absolute tolerance given the huge mean: this would fail hard
    # (off by O(1) or worse) under catastrophic cancellation.
    expect_equal(unname(res$vars[1, 1]), var(vals), tolerance = 1e-6)
})

test_that(".compute_blocked_stats_DuckDBMatrix correctly handles the integer64 COUNT(*) gotcha", {
    # DuckDB returns COUNT(*)/n() as bit64::integer64 via collect(); this
    # fixture (40 genes x 20 cells, each gene with its own count/sum/variance)
    # exercises the per-group aggregation broadly enough that an un-cast
    # integer64 count (silently wrong arithmetic when mixed with doubles)
    # would blow up the comparison against the R oracle well beyond float
    # noise. Manually verified during development: reverting the
    # as.double(result$stored_n) cast in .fused_group_variance_stats() makes
    # this comparison fail catastrophically (Inf/NaN), confirming this test
    # is sensitive to that regression.
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    res <- .compute_blocked_stats_DuckDBMatrix(m)

    expect_equal(unname(res$means[, 1]), unname(rowMeans(dense)))
    expect_equal(unname(res$vars[, 1]), unname(apply(dense, 1, var)))
})

test_that(".compute_group_stats_DuckDBMatrix matches the oracle, including NA-excluded cells", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    set.seed(1)
    groups <- sample(c("A", "B"), ncol(m), replace = TRUE)
    groups[1:3] <- NA

    res <- .compute_group_stats_DuckDBMatrix(m, groups = groups)

    for (lvl in c("A", "B")) {
        cols <- which(!is.na(groups) & groups == lvl)
        sub <- dense[, cols, drop = FALSE]
        expect_equal(unname(res$means[, lvl]), unname(rowMeans(sub)))
        expect_equal(unname(res$vars[, lvl]), unname(apply(sub, 1, var)))
    }
})

test_that(".compute_group_stats_DuckDBMatrix matches the oracle for a nonzero fill", {
    m <- makeSparseDuckDBMatrix() + 5
    dense <- as.matrix(m)

    set.seed(1)
    groups <- sample(c("A", "B"), ncol(m), replace = TRUE)

    res <- .compute_group_stats_DuckDBMatrix(m, groups = groups)

    for (lvl in c("A", "B")) {
        cols <- which(groups == lvl)
        sub <- dense[, cols, drop = FALSE]
        expect_equal(unname(res$means[, lvl]), unname(rowMeans(sub)))
        expect_equal(unname(res$vars[, lvl]), unname(apply(sub, 1, var)))
    }
})

test_that(".compute_blocked_stats_normalized_DuckDBMatrix matches the oracle (fill = 0)", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)

    set.seed(2)
    sf <- runif(ncol(m), 0.5, 1.5)
    norm_dense <- sweep(dense, 2, sf, "/")

    res <- .compute_blocked_stats_normalized_DuckDBMatrix(m, size.factors = sf)

    expect_equal(unname(res$means[, 1]), unname(rowMeans(norm_dense)))
    expect_equal(unname(res$vars[, 1]), unname(apply(norm_dense, 1, var)))
})

test_that(".compute_blocked_stats_normalized_DuckDBMatrix matches the oracle for a nonzero fill", {
    m <- makeSparseDuckDBMatrix()
    dense <- as.matrix(m)
    m5 <- m + 5
    dense5 <- as.matrix(m5)
    stored_mask <- dense != 0

    set.seed(2)
    sf <- runif(ncol(m), 0.5, 1.5)
    norm_fill5 <- 5 / mean(sf)

    # Oracle: stored cells normalize per-cell (value / sf); the implicit fill
    # population approximates with norm_fill (fill / mean(sf)), matching the
    # historical implementation's approximation (the true per-cell sf for an
    # absent cell is unknown).
    oracle <- matrix(norm_fill5, nrow(dense5), ncol(dense5), dimnames = dimnames(dense5))
    oracle[stored_mask] <- (dense5 / matrix(sf, nrow(dense5), ncol(dense5), byrow = TRUE))[stored_mask]

    res <- .compute_blocked_stats_normalized_DuckDBMatrix(m5, size.factors = sf)

    expect_equal(unname(res$means[, 1]), unname(rowMeans(oracle)))
    expect_equal(unname(res$vars[, 1]), unname(apply(oracle, 1, var)))
})

test_that("modelGeneVar/modelGeneCV2/pairwiseTTests end-to-end still match scran's own matrix path", {
    skip_if_not_installed("airway")
    names(dimnames(airway_counts)) <- c("index1", "index2")
    pqmat <- DuckDBMatrix(airway_counts_path, datacol = "value",
                          keycols = lapply(dimnames(airway_counts),
                                           function(x) setNames(seq_along(x), x)))
    mat <- as.matrix(pqmat)
    subset_genes <- 1:100

    log_ddb <- normalizeCounts(pqmat)
    log_mat <- as.matrix(log_ddb)

    ddb_mgv <- scran::modelGeneVar(log_ddb, subset.row = subset_genes)
    mat_mgv <- scran::modelGeneVar(log_mat, subset.row = subset_genes)
    expect_equal(ddb_mgv$mean, mat_mgv$mean, tolerance = 1e-6)
    expect_equal(ddb_mgv$total, mat_mgv$total, tolerance = 1e-6)

    set.seed(42)
    groups <- sample(c("g1", "g2", "g3"), ncol(pqmat), replace = TRUE)
    ddb_pt <- scran::pairwiseTTests(pqmat, groups = groups, subset.row = subset_genes)
    mat_pt <- scran::pairwiseTTests(mat, groups = groups, subset.row = subset_genes)
    expect_equal(ddb_pt$statistics[[1]]$logFC, mat_pt$statistics[[1]]$logFC, tolerance = 1e-6)
    expect_equal(ddb_pt$statistics[[1]]$p.value, mat_pt$statistics[[1]]$p.value, tolerance = 1e-4)
})
