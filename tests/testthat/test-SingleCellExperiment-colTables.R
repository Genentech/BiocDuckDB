# Checks for proper functioning of the colTable methods.
# library(BiocDuckDB); library(testthat); source("setup.R"); source("test-SingleCellExperiment-colTables.R")

library(SingleCellExperiment)

# Set up test data
set.seed(1000)
ncells <- 100L
nfeatures <- 200L
u <- matrix(rpois(nfeatures * ncells, 5), nrow=nfeatures, ncol=ncells)
v <- matrix(rnorm(nfeatures * ncells), nrow=nfeatures, ncol=ncells)
sce <- SingleCellExperiment(list(counts=u, logcounts=v))

# Mock up some nested sample tables (DataFrames with ncells rows)
t1 <- DataFrame(
    disease_id = paste0("MONDO:", sample(100000:999999, ncells, replace=TRUE)),
    disease_label = paste0("Disease_", 1:ncells),
    onset_age = sample(20:80, ncells, replace=TRUE)
)
t2 <- DataFrame(
    qc_metric = sample(c("n_genes", "n_counts", "pct_mito"), ncells, replace=TRUE),
    value = runif(ncells, 1000, 5000),
    passed = sample(c(TRUE, FALSE), ncells, replace=TRUE)
)

test_that("colTable getters/setters are functioning with character 'type'", {
    colTable(sce, "diseases") <- t1
    expect_identical(colTable(sce, "diseases"), t1)
    expect_identical(colTables(sce), SimpleList(diseases=t1))
    expect_identical(colTableNames(sce), "diseases")

    colTable(sce, "qc_metrics") <- t2
    expect_identical(colTable(sce, "qc_metrics"), t2)
    expect_identical(colTables(sce), SimpleList(diseases=t1, qc_metrics=t2))
    expect_identical(colTableNames(sce), c("diseases", "qc_metrics"))

    # Clearing values.
    colTable(sce, "diseases") <- NULL
    expect_identical(colTable(sce, "qc_metrics"), t2)
    expect_identical(colTables(sce), SimpleList(qc_metrics=t2))
    expect_identical(colTableNames(sce), "qc_metrics")

    # Checking for different errors.
    expect_error(colTable(sce, "diseases"), "invalid subscript")
    expect_error(colTable(sce, 2), "invalid subscript")
    expect_error(colTable(sce, "new") <- t1[1:10,], "number of rows")
    expect_error(colTable(sce, 1) <- "huh", "number of rows")
})

test_that("colTables getters/setters are functioning", {
    colTables(sce) <- list(diseases=t1, qc_metrics=t2)
    expect_identical(colTableNames(sce), c("diseases", "qc_metrics"))
    expect_identical(colTable(sce), t1)
    expect_identical(colTable(sce, "diseases"), t1)
    expect_identical(colTable(sce, 1), t1)
    expect_identical(colTable(sce, "qc_metrics"), t2)
    expect_identical(colTable(sce, 2), t2)

    # Clearing via empty List.
    alt <- sce
    colTables(alt) <- SimpleList()
    expect_identical(colTables(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    colTables(sce) <- SimpleList(patient_info=t1)
    expect_identical(SimpleList(patient_info=t1), colTables(sce))
    expect_identical(t1, colTable(sce,1))

    alt <- sce
    colTables(alt) <- NULL
    expect_identical(colTables(alt), setNames(SimpleList(), character(0)))

    # Checking for errors.
    expect_error(colTables(sce) <- list(t1, t2[1:10,]), "number of rows")
    expect_error(colTables(sce) <- list(t1[1:10,], t2[1:10,]), "number of rows")
})

test_that("getters/setters respond to dimnames", {
    named <- sce
    colnames(named) <- paste0("Cell", seq_len(ncol(named)))

    colTable(named, "diseases") <- t1
    colTable(named, "qc_metrics") <- t2
    expect_identical(rownames(colTable(named, 1)), colnames(named))
    expect_identical(rownames(colTable(named, 2)), colnames(named))
    expect_identical(rownames(colTable(named, 1, withDimnames=FALSE)), NULL)

    out <- colTables(named)
    expect_identical(rownames(out[[1]]), colnames(named))
    expect_identical(rownames(out[[2]]), colnames(named))
    out <- colTables(named, withDimnames=FALSE)
    expect_identical(rownames(out[[1]]), NULL)
    expect_identical(rownames(out[[2]]), NULL)

    # withDimnames works on the left hand side.
    rownames(colTable(named, "diseases", withDimnames=FALSE)) <- toupper(colnames(named))
    expect_identical(rownames(colTable(named, 1)), colnames(named))
    expect_identical(rownames(colTable(named, 1, withDimnames=FALSE)), toupper(colnames(named)))

    names(colTables(named, withDimnames=FALSE)) <- c("alpha", "bravo")
    expect_identical(rownames(colTable(named, 1, withDimnames=FALSE)), toupper(colnames(named)))

    # No warning when names are the same.
    t1.2 <- t1
    rownames(t1.2) <- colnames(named)
    colTable(named, 1) <- t1.2

    # withDimnames doesn't raise warnings on non-identity (different from reducedDims).
    t1.2 <- t1
    rownames(t1.2) <- toupper(colnames(named))
    expect_warning(colTable(named, "diseases") <- t1.2, "should be the same")
    expect_warning(colTables(named) <- list(diseases=t1.2), "should be the same")
    colTable(named, "diseases") <- t1
})

test_that("colTables getters/setters preserve mcols and metadata", {
    stuff <- List(diseases=t1, qc_metrics=t2)
    mcols(stuff)$A <- c("one", "two")
    metadata(stuff)$B <- "three"

    colTables(sce) <- stuff
    out <- colTables(sce)
    expect_identical(mcols(out), mcols(stuff))
    expect_identical(metadata(out), metadata(stuff))
})

test_that("colTable getters/setters work with numeric indices", {
    empty_sce <- SingleCellExperiment(list(counts=u))
    expect_error(colTable(empty_sce, 2), "invalid subscript")
    expect_error(colTable(empty_sce, "diseases"), "invalid subscript")

    expect_error(colTable(empty_sce, 1) <- t1, "out of bounds")
    expect_error(colTable(empty_sce, 2) <- t1, "out of bounds")

    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    colTables(empty_sce) <- list(t1, t2)
    expect_identical(colTable(empty_sce, 1), t1)
    expect_identical(colTable(empty_sce, 2), t2)

    t1_mod <- t1
    t1_mod$new_col <- seq_len(nrow(t1_mod))
    colTable(empty_sce, "diseases") <- t1_mod # t1 is still the first element.
    expect_identical(colTable(empty_sce, 1), t1)
    expect_identical(colTable(empty_sce, 2), t2)
    expect_identical(colTable(empty_sce, 3), t1_mod)

    colTable(empty_sce, 1) <- NULL # t2 becomes the first element now.
    expect_identical(colTable(empty_sce, 1), t2)
    expect_identical(colTable(empty_sce, 2), colTable(empty_sce, "diseases"))

    colTable(empty_sce, 1) <- NULL # 't1_mod' becomes the first element.
    expect_identical(colTable(empty_sce, 1), t1_mod)
    expect_identical(colTableNames(empty_sce), "diseases")
    colTable(empty_sce, 1) <- t2 # t2 now overwrites the first element.
    expect_identical(colTable(empty_sce, 1), t2)
    expect_identical(colTableNames(empty_sce), "diseases")

    expect_error(colTable(empty_sce, 5) <- t1, "out of bounds")
})

test_that("colTableNames getters/setters work correctly", {
    colTables(sce) <- list(t1, t2)
    expect_true(all(c("unnamed1", "unnamed2") %in% colTableNames(sce) | 
                    c("1", "2") %in% colTableNames(sce)))

    colTables(sce) <- list(diseases=t1, qc_metrics=t2)
    expect_identical(colTableNames(sce), c("diseases", "qc_metrics"))

    # Directly setting.
    colTableNames(sce) <- c("A", "B")
    expect_identical(colTableNames(sce), c("A", "B"))

    # When wiped.
    colTables(sce) <- NULL
    expect_identical(colTableNames(sce), character(0))

    # Setting names with wrong length should error
    colTables(sce) <- list(diseases=t1, qc_metrics=t2)
    expect_error(colTableNames(sce) <- c("A", "B", "C"), "more column names")
})

test_that("colTables work after subsetting SCE by columns", {
    colnames(sce) <- paste0("Cell", seq_len(ncol(sce)))
    colTable(sce, "diseases") <- t1
    colTable(sce, "qc_metrics") <- t2

    # Subset to first 50 cells
    sce_sub <- sce[, 1:50]

    # Check tables are subset correctly
    expect_identical(nrow(colTable(sce_sub)), 50L)
    expect_identical(nrow(colTable(sce_sub, "diseases")), 50L)
    expect_identical(nrow(colTable(sce_sub, "qc_metrics")), 50L)
    expect_identical(rownames(colTable(sce_sub)), colnames(sce_sub))
    expect_identical(rownames(colTable(sce_sub, "diseases")), colnames(sce_sub))
    expect_identical(rownames(colTable(sce_sub, "qc_metrics")), colnames(sce_sub))

    # Check values are correct
    expect_identical(colTable(sce_sub, withDimnames=FALSE), t1[1:50, ])
    expect_identical(colTable(sce_sub, "diseases", withDimnames=FALSE), t1[1:50, ])
    expect_identical(colTable(sce_sub, "qc_metrics", withDimnames=FALSE), t2[1:50, ])
})

test_that("colTables are independent from reducedDims", {
    # Add both reducedDims and colTables
    d1 <- matrix(rnorm(ncells * 4), ncol=4)  # cells × components
    reducedDim(sce, "PCA") <- d1
    colTable(sce, "diseases") <- t1

    # Both should coexist
    expect_identical(reducedDimNames(sce), "PCA")
    expect_identical(colTableNames(sce), "diseases")
    expect_true(is.matrix(reducedDim(sce, "PCA")))
    expect_s4_class(colTable(sce, "diseases"), "DFrame")
})

test_that("colTables handle data.frame conversion", {
    # Test with regular data.frame (should convert to DataFrame)
    df1 <- data.frame(
        col1 = paste0("val", 1:ncells),
        col2 = rnorm(ncells)
    )
    colTable(sce, "from_df") <- df1
    expect_s4_class(colTable(sce, "from_df"), "DFrame")
    expect_identical(nrow(colTable(sce, "from_df")), ncells)
})

test_that("colTables storage in int_colData", {
    colTable(sce, "diseases") <- t1

    # Check it's stored in the right place
    expect_true("colTables" %in% names(int_colData(sce)))

    # Check structure
    internal <- int_colData(sce)[["colTables"]]
    expect_s4_class(internal, "DFrame")
    expect_identical(nrow(internal), ncol(sce))
})
