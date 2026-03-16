# Checks for proper functioning of the rowLoading methods.
# library(BiocDuckDB); library(testthat); source("setup.R"); source("test-SingleCellExperiment-loadings.R")

library(SingleCellExperiment)

# Set up test data
set.seed(1000)
ncells <- 100
nfeatures <- 200
u <- matrix(rpois(nfeatures * ncells, 5), nrow=nfeatures, ncol=ncells)
v <- matrix(rnorm(nfeatures * ncells), nrow=nfeatures, ncol=ncells)
sce <- SingleCellExperiment(list(counts=u, logcounts=v))

# Mock up some loading matrices (features × components)
l1 <- matrix(rnorm(nfeatures * 5), nrow=nfeatures, ncol=5)  # PCA loadings: 200 genes × 5 PCs
l2 <- matrix(rnorm(nfeatures * 20), nrow=nfeatures, ncol=20)  # Factor loadings: 200 genes × 20 factors
colnames(l1) <- paste0("PC", 1:5)
colnames(l2) <- paste0("Factor", 1:20)

test_that("rowLoading getters/setters are functioning with character 'type'", {
    rowLoading(sce, "PCA") <- l1
    expect_identical(rowLoading(sce, "PCA"), l1)
    expect_identical(rowLoadings(sce), SimpleList(PCA=l1))
    expect_identical(rowLoadingNames(sce), "PCA")

    rowLoading(sce, "Factors") <- l2
    expect_identical(rowLoading(sce, "Factors"), l2)
    expect_identical(rowLoadings(sce), SimpleList(PCA=l1, Factors=l2))
    expect_identical(rowLoadingNames(sce), c("PCA", "Factors"))

    # Clearing values.
    rowLoading(sce, "PCA") <- NULL
    expect_identical(rowLoading(sce, "Factors"), l2)
    expect_identical(rowLoadings(sce), SimpleList(Factors=l2))
    expect_identical(rowLoadingNames(sce), "Factors")

    # Checking for different errors.
    expect_error(rowLoading(sce, "PCA"), "invalid subscript")
    expect_error(rowLoading(sce, 2), "invalid subscript")
    expect_error(rowLoading(sce, "ICA") <- l1[1:10,], "number of rows")
    expect_error(rowLoading(sce, 1) <- "huh", "number of rows")
})

test_that("rowLoadings getters/setters are functioning", {
    rowLoadings(sce) <- list(PCA=l1, Factors=l2)
    expect_identical(rowLoadingNames(sce), c("PCA", "Factors"))
    expect_identical(rowLoading(sce), l1)
    expect_identical(rowLoading(sce, "PCA"), l1)
    expect_identical(rowLoading(sce, 1), l1)
    expect_identical(rowLoading(sce, "Factors"), l2)
    expect_identical(rowLoading(sce, 2), l2)

    # Clearing via empty List.
    alt <- sce
    rowLoadings(alt) <- SimpleList()
    expect_identical(rowLoadings(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    rowLoadings(sce) <- SimpleList(ICA=l1)
    expect_identical(SimpleList(ICA=l1), rowLoadings(sce))
    expect_identical(l1, rowLoading(sce,1))

    alt <- sce
    rowLoadings(alt) <- NULL
    expect_identical(rowLoadings(alt), setNames(SimpleList(), character(0)))

    # Checking for errors.
    expect_error(rowLoadings(sce) <- list(l1, l2[1:10,]), "number of rows")
    expect_error(rowLoadings(sce) <- list(l1[1:10,], l2[1:10,]), "number of rows")
})

test_that("getters/setters respond to dimnames", {
    named <- sce
    rownames(named) <- paste0("Gene", seq_len(nrow(named)))

    rowLoading(named, "PCA") <- l1
    rowLoading(named, "Factors") <- l2
    expect_identical(rownames(rowLoading(named, 1)), rownames(named))
    expect_identical(rownames(rowLoading(named, 2)), rownames(named))
    expect_identical(rownames(rowLoading(named, 1, withDimnames=FALSE)), NULL)

    out <- rowLoadings(named)
    expect_identical(rownames(out[[1]]), rownames(named))
    expect_identical(rownames(out[[2]]), rownames(named))
    out <- rowLoadings(named, withDimnames=FALSE)
    expect_identical(rownames(out[[1]]), NULL)
    expect_identical(rownames(out[[2]]), NULL)

    # withDimnames works on the left hand side.
    rownames(rowLoading(named, "PCA", withDimnames=FALSE)) <- toupper(rownames(named))
    expect_identical(rownames(rowLoading(named, 1)), rownames(named))
    expect_identical(rownames(rowLoading(named, 1, withDimnames=FALSE)), toupper(rownames(named)))

    names(rowLoadings(named, withDimnames=FALSE)) <- c("alpha", "bravo")
    expect_identical(rownames(rowLoading(named, 1, withDimnames=FALSE)), toupper(rownames(named)))

    # No warning when names are the same.
    l1.2 <- l1
    rownames(l1.2) <- rownames(named)
    rowLoading(named, 1) <- l1.2

    # withDimnames doesn't raise warnings on non-identity (different from reducedDims).
    l1.2 <- l1
    rownames(l1.2) <- toupper(rownames(named))
    expect_warning(rowLoading(named, "PCA") <- l1.2, "should be the same")
    expect_warning(rowLoadings(named) <- list(PCA=l1.2), "should be the same")
    rowLoading(named, "PCA") <- l1
})

test_that("rowLoadings getters/setters preserve mcols and metadata", {
    stuff <- List(PCA=l1, Factors=l2)
    mcols(stuff)$A <- c("one", "two")
    metadata(stuff)$B <- "three"

    rowLoadings(sce) <- stuff
    out <- rowLoadings(sce)
    expect_identical(mcols(out), mcols(stuff))
    expect_identical(metadata(out), metadata(stuff))
})

test_that("rowLoading getters/setters work with numeric indices", {
    empty_sce <- SingleCellExperiment(list(counts=u))
    expect_error(rowLoading(empty_sce, 2), "invalid subscript") 
    expect_error(rowLoading(empty_sce, "PCA"), "invalid subscript") 

    expect_error(rowLoading(empty_sce, 1) <- l1, "out of bounds")
    expect_error(rowLoading(empty_sce, 2) <- l1, "out of bounds")

    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    rowLoadings(empty_sce) <- list(l1, l2)
    expect_identical(rowLoading(empty_sce, 1), l1)
    expect_identical(rowLoading(empty_sce, 2), l2)

    mult <- l1 * 5
    rowLoading(empty_sce, "PCA") <- mult # l1 is still the first element.
    expect_identical(rowLoading(empty_sce, 1), l1)
    expect_identical(rowLoading(empty_sce, 2), l2)
    expect_identical(rowLoading(empty_sce, 3), mult)

    rowLoading(empty_sce, 1) <- NULL # l2 becomes the first element now.
    expect_identical(rowLoading(empty_sce, 1), l2)
    expect_identical(rowLoading(empty_sce, 2), rowLoading(empty_sce, "PCA"))

    rowLoading(empty_sce, 1) <- NULL # 'mult' becomes the first element.
    expect_identical(rowLoading(empty_sce, 1), mult)
    expect_identical(rowLoadingNames(empty_sce), "PCA")
    rowLoading(empty_sce, 1) <- l2 # l2 now overwrites the first element.
    expect_identical(rowLoading(empty_sce, 1), l2)
    expect_identical(rowLoadingNames(empty_sce), "PCA")

    expect_error(rowLoading(empty_sce, 5) <- l1, "out of bounds")
})

test_that("rowLoadingNames getters/setters work correctly", {
    rowLoadings(sce) <- list(l1, l2)
    expect_true(all(c("unnamed1", "unnamed2") %in% rowLoadingNames(sce) | 
                    c("1", "2") %in% rowLoadingNames(sce)))

    rowLoadings(sce) <- list(PCA=l1, Factors=l2)
    expect_identical(rowLoadingNames(sce), c("PCA", "Factors"))

    # Directly setting.
    rowLoadingNames(sce) <- c("A", "B")
    expect_identical(rowLoadingNames(sce), c("A", "B"))

    # When wiped.
    rowLoadings(sce) <- NULL
    expect_identical(rowLoadingNames(sce), character(0))

    # Setting names with wrong length should error
    rowLoadings(sce) <- list(PCA=l1, Factors=l2)
    expect_error(rowLoadingNames(sce) <- c("A", "B", "C"), "more column names")
})

test_that("rowLoadings work after subsetting SCE by rows", {
    rownames(sce) <- paste0("Gene", seq_len(nrow(sce)))
    rowLoading(sce, "PCA") <- l1
    rowLoading(sce, "Factors") <- l2

    # Subset to first 100 genes
    sce_sub <- sce[1:100, ]

    # Check loadings are subset correctly
    expect_identical(nrow(rowLoading(sce_sub)), 100L)
    expect_identical(nrow(rowLoading(sce_sub, "PCA")), 100L)
    expect_identical(nrow(rowLoading(sce_sub, "Factors")), 100L)
    expect_identical(rownames(rowLoading(sce_sub)), rownames(sce_sub))
    expect_identical(rownames(rowLoading(sce_sub, "PCA")), rownames(sce_sub))
    expect_identical(rownames(rowLoading(sce_sub, "Factors")), rownames(sce_sub))

    # Check values are correct
    expect_identical(rowLoading(sce_sub, withDimnames=FALSE), l1[1:100, ])
    expect_identical(rowLoading(sce_sub, "PCA", withDimnames=FALSE), l1[1:100, ])
    expect_identical(rowLoading(sce_sub, "Factors", withDimnames=FALSE), l2[1:100, ])
})

test_that("rowLoadings are independent from reducedDims", {
    # Add both reducedDims and rowLoadings
    d1 <- matrix(rnorm(ncells * 4), ncol=4)  # cells × components
    reducedDim(sce, "PCA_embedding") <- d1
    rowLoading(sce, "PCA_loading") <- l1

    # Both should coexist
    expect_identical(reducedDimNames(sce), "PCA_embedding")
    expect_identical(rowLoadingNames(sce), "PCA_loading")
    expect_identical(nrow(reducedDim(sce, "PCA_embedding")), ncol(sce))
    expect_identical(nrow(rowLoading(sce, "PCA_loading")), nrow(sce))
})

test_that("rowLoadings handle different matrix types", {
    # Test with sparse matrix
    sparse_l1 <- Matrix::Matrix(l1, sparse=TRUE)
    rowLoading(sce, "sparse") <- sparse_l1
    expect_s4_class(rowLoading(sce, "sparse"), "Matrix")

    # Test with regular matrix
    rowLoading(sce, "dense") <- l1
    expect_true(is.matrix(rowLoading(sce, "dense")))
})

test_that("rowLoadings storage in int_elementMetadata", {
    rowLoading(sce, "PCA") <- l1

    # Check it's stored in the right place
    expect_true("rowLoadings" %in% names(int_elementMetadata(sce)))

    # Check structure
    internal <- int_elementMetadata(sce)[["rowLoadings"]]
    expect_s4_class(internal, "DFrame")
    expect_identical(nrow(internal), nrow(sce))
})
