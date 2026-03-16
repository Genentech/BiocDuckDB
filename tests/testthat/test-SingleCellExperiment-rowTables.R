# Checks for proper functioning of the rowTable methods.
# library(BiocDuckDB); library(testthat); source("setup.R"); source("test-SingleCellExperiment-rowTables.R")

library(SingleCellExperiment)

# Set up test data
set.seed(1000)
ncells <- 100L
nfeatures <- 200L
u <- matrix(rpois(nfeatures * ncells, 5), nrow=nfeatures, ncol=ncells)
v <- matrix(rnorm(nfeatures * ncells), nrow=nfeatures, ncol=ncells)
sce <- SingleCellExperiment(list(counts=u, logcounts=v))

# Mock up some nested feature tables (DataFrames with nfeatures rows)
t1 <- DataFrame(
    isoform_id = paste0("ISO", 1:nfeatures),
    length = sample(500:5000, nfeatures, replace=TRUE),
    is_canonical = sample(c(TRUE, FALSE), nfeatures, replace=TRUE)
)
t2 <- DataFrame(
    go_term = paste0("GO:", sample(1000:9999, nfeatures, replace=TRUE)),
    ontology = sample(c("BP", "MF", "CC"), nfeatures, replace=TRUE),
    evidence = sample(c("IDA", "IEA", "ISS"), nfeatures, replace=TRUE)
)

test_that("rowTable getters/setters are functioning with character 'type'", {
    rowTable(sce, "isoforms") <- t1
    expect_identical(rowTable(sce, "isoforms"), t1)
    expect_identical(rowTables(sce), SimpleList(isoforms=t1))
    expect_identical(rowTableNames(sce), "isoforms")

    rowTable(sce, "go_terms") <- t2
    expect_identical(rowTable(sce, "go_terms"), t2)
    expect_identical(rowTables(sce), SimpleList(isoforms=t1, go_terms=t2))
    expect_identical(rowTableNames(sce), c("isoforms", "go_terms"))

    # Clearing values.
    rowTable(sce, "isoforms") <- NULL
    expect_identical(rowTable(sce, "go_terms"), t2)
    expect_identical(rowTables(sce), SimpleList(go_terms=t2))
    expect_identical(rowTableNames(sce), "go_terms")

    # Checking for different errors.
    expect_error(rowTable(sce, "isoforms"), "invalid subscript")
    expect_error(rowTable(sce, 2), "invalid subscript")
    expect_error(rowTable(sce, "new") <- t1[1:10,], "number of rows")
    expect_error(rowTable(sce, 1) <- "huh", "number of rows")
})

test_that("rowTables getters/setters are functioning", {
    rowTables(sce) <- list(isoforms=t1, go_terms=t2)
    expect_identical(rowTableNames(sce), c("isoforms", "go_terms"))
    expect_identical(rowTable(sce), t1)
    expect_identical(rowTable(sce, "isoforms"), t1)
    expect_identical(rowTable(sce, 1), t1)
    expect_identical(rowTable(sce, "go_terms"), t2)
    expect_identical(rowTable(sce, 2), t2)

    # Clearing via empty List.
    alt <- sce
    rowTables(alt) <- SimpleList()
    expect_identical(rowTables(alt), setNames(SimpleList(), character(0)))

    # Clearing via NULL.
    rowTables(sce) <- SimpleList(de_stats=t1)
    expect_identical(SimpleList(de_stats=t1), rowTables(sce))
    expect_identical(t1, rowTable(sce,1))

    alt <- sce
    rowTables(alt) <- NULL
    expect_identical(rowTables(alt), setNames(SimpleList(), character(0)))

    # Checking for errors.
    expect_error(rowTables(sce) <- list(t1, t2[1:10,]), "number of rows")
    expect_error(rowTables(sce) <- list(t1[1:10,], t2[1:10,]), "number of rows")
})

test_that("getters/setters respond to dimnames", {
    named <- sce
    rownames(named) <- paste0("Gene", seq_len(nrow(named)))

    rowTable(named, "isoforms") <- t1
    rowTable(named, "go_terms") <- t2
    expect_identical(rownames(rowTable(named, 1)), rownames(named))
    expect_identical(rownames(rowTable(named, 2)), rownames(named))
    expect_identical(rownames(rowTable(named, 1, withDimnames=FALSE)), NULL)

    out <- rowTables(named)
    expect_identical(rownames(out[[1]]), rownames(named))
    expect_identical(rownames(out[[2]]), rownames(named))
    out <- rowTables(named, withDimnames=FALSE)
    expect_identical(rownames(out[[1]]), NULL)
    expect_identical(rownames(out[[2]]), NULL)

    # withDimnames works on the left hand side.
    rownames(rowTable(named, "isoforms", withDimnames=FALSE)) <- toupper(rownames(named))
    expect_identical(rownames(rowTable(named, 1)), rownames(named))
    expect_identical(rownames(rowTable(named, 1, withDimnames=FALSE)), toupper(rownames(named)))

    names(rowTables(named, withDimnames=FALSE)) <- c("alpha", "bravo")
    expect_identical(rownames(rowTable(named, 1, withDimnames=FALSE)), toupper(rownames(named)))

    # No warning when names are the same.
    t1.2 <- t1
    rownames(t1.2) <- rownames(named)
    rowTable(named, 1) <- t1.2

    # withDimnames doesn't raise warnings on non-identity (different from reducedDims).
    t1.2 <- t1
    rownames(t1.2) <- toupper(rownames(named))
    expect_warning(rowTable(named, "isoforms") <- t1.2, "should be the same")
    expect_warning(rowTables(named) <- list(isoforms=t1.2), "should be the same")
    rowTable(named, "isoforms") <- t1
})

test_that("rowTables getters/setters preserve mcols and metadata", {
    stuff <- List(isoforms=t1, go_terms=t2)
    mcols(stuff)$A <- c("one", "two")
    metadata(stuff)$B <- "three"

    rowTables(sce) <- stuff
    out <- rowTables(sce)
    expect_identical(mcols(out), mcols(stuff))
    expect_identical(metadata(out), metadata(stuff))
})

test_that("rowTable getters/setters work with numeric indices", {
    empty_sce <- SingleCellExperiment(list(counts=u))
    expect_error(rowTable(empty_sce, 2), "invalid subscript")
    expect_error(rowTable(empty_sce, "isoforms"), "invalid subscript")

    expect_error(rowTable(empty_sce, 1) <- t1, "out of bounds")
    expect_error(rowTable(empty_sce, 2) <- t1, "out of bounds")

    # This gets a bit confusing as the order changes when earlier elements are wiped out.
    rowTables(empty_sce) <- list(t1, t2)
    expect_identical(rowTable(empty_sce, 1), t1)
    expect_identical(rowTable(empty_sce, 2), t2)

    t1_mod <- t1
    t1_mod$new_col <- seq_len(nrow(t1_mod))
    rowTable(empty_sce, "isoforms") <- t1_mod # t1 is still the first element.
    expect_identical(rowTable(empty_sce, 1), t1)
    expect_identical(rowTable(empty_sce, 2), t2)
    expect_identical(rowTable(empty_sce, 3), t1_mod)

    rowTable(empty_sce, 1) <- NULL # t2 becomes the first element now.
    expect_identical(rowTable(empty_sce, 1), t2)
    expect_identical(rowTable(empty_sce, 2), rowTable(empty_sce, "isoforms"))

    rowTable(empty_sce, 1) <- NULL # 't1_mod' becomes the first element.
    expect_identical(rowTable(empty_sce, 1), t1_mod)
    expect_identical(rowTableNames(empty_sce), "isoforms")
    rowTable(empty_sce, 1) <- t2 # t2 now overwrites the first element.
    expect_identical(rowTable(empty_sce, 1), t2)
    expect_identical(rowTableNames(empty_sce), "isoforms")

    expect_error(rowTable(empty_sce, 5) <- t1, "out of bounds")
})

test_that("rowTableNames getters/setters work correctly", {
    rowTables(sce) <- list(t1, t2)
    expect_true(all(c("unnamed1", "unnamed2") %in% rowTableNames(sce) | 
                    c("1", "2") %in% rowTableNames(sce)))

    rowTables(sce) <- list(isoforms=t1, go_terms=t2)
    expect_identical(rowTableNames(sce), c("isoforms", "go_terms"))

    # Directly setting.
    rowTableNames(sce) <- c("A", "B")
    expect_identical(rowTableNames(sce), c("A", "B"))

    # When wiped.
    rowTables(sce) <- NULL
    expect_identical(rowTableNames(sce), character(0))

    # Setting names with wrong length should error
    rowTables(sce) <- list(isoforms=t1, go_terms=t2)
    expect_error(rowTableNames(sce) <- c("A", "B", "C"), "more column names")
})

test_that("rowTables work after subsetting SCE by rows", {
    rownames(sce) <- paste0("Gene", seq_len(nrow(sce)))
    rowTable(sce, "isoforms") <- t1
    rowTable(sce, "go_terms") <- t2

    # Subset to first 100 genes
    sce_sub <- sce[1:100, ]

    # Check tables are subset correctly
    expect_identical(nrow(rowTable(sce_sub)), 100L)
    expect_identical(nrow(rowTable(sce_sub, "isoforms")), 100L)
    expect_identical(nrow(rowTable(sce_sub, "go_terms")), 100L)
    expect_identical(rownames(rowTable(sce_sub)), rownames(sce_sub))
    expect_identical(rownames(rowTable(sce_sub, "isoforms")), rownames(sce_sub))
    expect_identical(rownames(rowTable(sce_sub, "go_terms")), rownames(sce_sub))

    # Check values are correct
    expect_identical(rowTable(sce_sub, withDimnames=FALSE), t1[1:100, ])
    expect_identical(rowTable(sce_sub, "isoforms", withDimnames=FALSE), t1[1:100, ])
    expect_identical(rowTable(sce_sub, "go_terms", withDimnames=FALSE), t2[1:100, ])
})

test_that("rowTables are independent from rowLoadings", {
    # Add both rowLoadings and rowTables
    l1 <- matrix(rnorm(nfeatures * 5), nrow=nfeatures, ncol=5)  # PCA loadings
    rowLoading(sce, "PCA") <- l1
    rowTable(sce, "isoforms") <- t1

    # Both should coexist
    expect_identical(rowLoadingNames(sce), "PCA")
    expect_identical(rowTableNames(sce), "isoforms")
    expect_true(is.matrix(rowLoading(sce, "PCA")))
    expect_s4_class(rowTable(sce, "isoforms"), "DFrame")
})

test_that("rowTables handle data.frame conversion", {
    # Test with regular data.frame (should convert to DataFrame)
    df1 <- data.frame(
        col1 = paste0("val", 1:nfeatures),
        col2 = rnorm(nfeatures)
    )
    rowTable(sce, "from_df") <- df1
    expect_s4_class(rowTable(sce, "from_df"), "DFrame")
    expect_identical(nrow(rowTable(sce, "from_df")), nfeatures)
})

test_that("rowTables storage in int_elementMetadata", {
    rowTable(sce, "isoforms") <- t1

    # Check it's stored in the right place
    expect_true("rowTables" %in% names(int_elementMetadata(sce)))

    # Check structure
    internal <- int_elementMetadata(sce)[["rowTables"]]
    expect_s4_class(internal, "DFrame")
    expect_identical(nrow(internal), nrow(sce))
})
