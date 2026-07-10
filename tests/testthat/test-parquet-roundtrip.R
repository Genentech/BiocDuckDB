# Comprehensive tests for writeParquet/readParquet round-trip functionality
# library(BiocDuckDB); library(testthat); source("setup.R"); source("test-parquet-roundtrip.R")

library(SingleCellExperiment)
library(MultiAssayExperiment)

# ==============================================================================
# SummarizedExperiment Tests
# ==============================================================================

test_that("SummarizedExperiment basic round-trip works", {
    set.seed(100)
    ncells <- 50L
    ngenes <- 100L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    logcounts <- log1p(counts)

    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    se <- SummarizedExperiment(
        assays = list(counts = counts, logcounts = logcounts),
        rowData = DataFrame(
            gene_type = sample(c("protein_coding", "lncRNA"), ngenes, replace = TRUE),
            chromosome = sample(paste0("chr", 1:22), ngenes, replace = TRUE)
        ),
        colData = DataFrame(
            sample_id = paste0("S", 1:ncells),
            batch = sample(c("A", "B", "C"), ncells, replace = TRUE),
            n_genes = rpois(ncells, 500)
        ),
        metadata = list(
            version = "1.0",
            date = "2026-03-24",
            description = "Test dataset"
        )
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(se, tmpdir)
    se2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(dim(se2), dim(se))
    expect_identical(assayNames(se2), assayNames(se))

    # Check assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_expected <- assay(se, "counts")
    logcounts_expected <- assay(se, "logcounts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    names(dimnames(logcounts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(se2, "counts"), counts_expected)
    checkDuckDBMatrix(assay(se2, "logcounts"), logcounts_expected)

    # Check row/col data (convert expected to data.frame)
    checkDuckDBDataFrame(rowData(se2), as.data.frame(rowData(se)))
    checkDuckDBDataFrame(colData(se2), as.data.frame(colData(se)))

    # Check metadata
    expect_identical(metadata(se2), metadata(se))

    unlink(tmpdir, recursive = TRUE)
})

test_that("RangedSummarizedExperiment with GenomicRanges works", {
    set.seed(101)
    ncells <- 30L
    ngenes <- 50L

    counts <- matrix(rpois(ngenes * ncells, 10), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    gr <- GRanges(
        seqnames = sample(paste0("chr", 1:5), ngenes, replace = TRUE),
        ranges = IRanges(
            start = sample(1:1000000, ngenes),
            width = sample(1000:5000, ngenes)
        ),
        strand = sample(c("+", "-"), ngenes, replace = TRUE),
        gene_id = paste0("ENSG", sprintf("%011d", 1:ngenes)),
        gene_name = paste0("Gene", seq_len(ngenes))
    )
    names(gr) <- rownames(counts)

    rse <- SummarizedExperiment(
        assays = list(counts = counts),
        rowRanges = gr,
        colData = DataFrame(cell_type = sample(c("A", "B"), ncells, replace = TRUE))
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(rse, tmpdir)
    rse2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(rse2, "RangedSummarizedExperiment")
    expect_identical(dim(rse2), dim(rse))

    # Check assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_expected <- assay(rse, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(rse2, "counts"), counts_expected)

    # Check rowRanges
    checkDuckDBGRanges(rowRanges(rse2), gr)

    # Check colData
    checkDuckDBDataFrame(colData(rse2), as.data.frame(colData(rse)))

    unlink(tmpdir, recursive = TRUE)
})

# ==============================================================================
# SingleCellExperiment Tests
# ==============================================================================

test_that("SingleCellExperiment with reducedDims round-trip works", {
    set.seed(102)
    ncells <- 40L
    ngenes <- 80L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    pca_mat <- matrix(rnorm(ncells * 10), ncol = 10)
    rownames(pca_mat) <- colnames(counts)

    umap_mat <- matrix(rnorm(ncells * 2), ncol = 2)
    rownames(umap_mat) <- colnames(counts)

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(cluster = sample(1:3, ncells, replace = TRUE)),
        reducedDims = list(PCA = pca_mat, UMAP = umap_mat)
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)
    sce2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_identical(dim(sce2), dim(sce))

    # Check assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_expected <- assay(sce, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(sce2, "counts"), counts_expected)

    # Check colData
    checkDuckDBDataFrame(colData(sce2), as.data.frame(colData(sce)))

    # Check reducedDims
    expect_identical(reducedDimNames(sce2), c("PCA", "UMAP"))
    checkDuckDBEmbeddings(reducedDim(sce2, "PCA"), pca_mat)
    checkDuckDBEmbeddings(reducedDim(sce2, "UMAP"), umap_mat)

    unlink(tmpdir, recursive = TRUE)
})

test_that("SingleCellExperiment with altExps (nested_experiment) works", {
    set.seed(103)
    ncells <- 30L
    ngenes_rna <- 50L
    ngenes_adt <- 20L

    # Main experiment (RNA)
    counts_rna <- matrix(rpois(ngenes_rna * ncells, 5), nrow = ngenes_rna, ncol = ncells)
    rownames(counts_rna) <- paste0("Gene", seq_len(ngenes_rna))
    colnames(counts_rna) <- paste0("Cell", seq_len(ncells))

    # ADT experiment
    counts_adt <- matrix(rpois(ngenes_adt * ncells, 100), nrow = ngenes_adt, ncol = ncells)
    rownames(counts_adt) <- paste0("Protein", seq_len(ngenes_adt))
    colnames(counts_adt) <- paste0("Cell", seq_len(ncells))

    adt_sce <- SingleCellExperiment(
        assays = list(counts = counts_adt),
        rowData = DataFrame(protein_id = paste0("P", 1:ngenes_adt))
    )

    sce <- SingleCellExperiment(
        assays = list(counts = counts_rna),
        colData = DataFrame(cell_type = sample(c("A", "B"), ncells, replace = TRUE)),
        altExps = list(ADT = adt_sce)
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)
    sce2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_identical(dim(sce2), dim(sce))

    # Check main assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_rna_expected <- assay(sce, "counts")
    names(dimnames(counts_rna_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(sce2, "counts"), counts_rna_expected)

    # Check colData
    checkDuckDBDataFrame(colData(sce2), as.data.frame(colData(sce)))

    # Check altExps
    expect_identical(altExpNames(sce2), "ADT")
    adt2 <- altExp(sce2, "ADT")
    expect_s4_class(adt2, "SingleCellExperiment")
    counts_adt_expected <- assay(altExp(sce, "ADT"), "counts")
    names(dimnames(counts_adt_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(adt2, "counts"), counts_adt_expected)
    checkDuckDBDataFrame(rowData(adt2), as.data.frame(rowData(adt_sce)))

    unlink(tmpdir, recursive = TRUE)
})

test_that("SingleCellExperiment with colTables works", {
    set.seed(104)
    ncells <- 25L
    ngenes <- 40L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(sample_id = paste0("S", 1:ncells))
    )

    # Add nested sample table with rownames
    disease_df <- DataFrame(
        disease_id = paste0("MONDO:", sample(100000:999999, ncells, replace = TRUE)),
        disease_label = paste0("Disease_", 1:ncells)
    )
    rownames(disease_df) <- colnames(counts)
    colTable(sce, "diseases") <- disease_df

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)
    sce2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_identical(dim(sce2), dim(sce))

    # Check assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_expected <- assay(sce, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(sce2, "counts"), counts_expected)

    # Check colData
    checkDuckDBDataFrame(colData(sce2), as.data.frame(colData(sce)))

    # Check colTables
    expect_identical(colTableNames(sce2), "diseases")
    checkDuckDBDataFrame(colTable(sce2, "diseases"), as.data.frame(disease_df))

    unlink(tmpdir, recursive = TRUE)
})

test_that("SingleCellExperiment with colPairs (graph_edges) works", {
    set.seed(105)
    ncells <- 20L
    ngenes <- 30L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        colData = DataFrame(cluster = sample(1:3, ncells, replace = TRUE))
    )

    # Add KNN graph
    knn_edges <- sample(1:ncells, ncells * 3, replace = TRUE)
    knn_hits <- SelfHits(
        from = rep(1:ncells, each = 3),
        to = knn_edges,
        nnode = ncells
    )
    colPair(sce, "knn") <- knn_hits

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)
    sce2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(sce2, "SingleCellExperiment")
    expect_identical(dim(sce2), dim(sce))

    # Check assays (add names(dimnames) since DuckDBMatrix adds them)
    counts_expected <- assay(sce, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(sce2, "counts"), counts_expected)

    # Check colData
    checkDuckDBDataFrame(colData(sce2), as.data.frame(colData(sce)))

    # Check colPairs
    expect_identical(names(colPairs(sce2)), "knn")
    checkDuckDBSelfHits(colPair(sce2, "knn"), colPair(sce, "knn"))

    unlink(tmpdir, recursive = TRUE)
})

# ==============================================================================
# MultiAssayExperiment Tests
# ==============================================================================

test_that("MultiAssayExperiment with nested experiments works", {
    set.seed(106)
    ncells <- 20L
    ngenes_rna <- 30L
    ngenes_protein <- 40L

    # Two separate experiments with proper structure
    counts_rna <- matrix(rpois(ngenes_rna * ncells, 5), nrow = ngenes_rna, ncol = ncells)
    rownames(counts_rna) <- paste0("Gene", seq_len(ngenes_rna))
    colnames(counts_rna) <- paste0("Cell", seq_len(ncells))

    counts_protein <- matrix(rpois(ngenes_protein * ncells, 5), nrow = ngenes_protein, ncol = ncells)
    rownames(counts_protein) <- paste0("Protein", seq_len(ngenes_protein))
    colnames(counts_protein) <- paste0("Cell", seq_len(ncells))

    sce1 <- SingleCellExperiment(assays = list(counts = counts_rna))
    sce2 <- SingleCellExperiment(assays = list(counts = counts_protein))

    # Create MAE with proper sample mapping
    mae <- MultiAssayExperiment(
        experiments = ExperimentList(RNA = sce1, Protein = sce2)
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(mae, tmpdir)
    mae2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(mae2, "MultiAssayExperiment")
    expect_identical(names(experiments(mae2)), c("RNA", "Protein"))
    expect_identical(length(experiments(mae2)), 2L)

    # Check each experiment
    expect_s4_class(mae2[["RNA"]], "SingleCellExperiment")
    expect_s4_class(mae2[["Protein"]], "SingleCellExperiment")

    # Check RNA experiment
    rna_expected <- assay(mae[["RNA"]], "counts")
    names(dimnames(rna_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(mae2[["RNA"]], "counts"), rna_expected)

    # Check Protein experiment
    protein_expected <- assay(mae[["Protein"]], "counts")
    names(dimnames(protein_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(mae2[["Protein"]], "counts"), protein_expected)

    unlink(tmpdir, recursive = TRUE)
})

test_that("MultiAssayExperiment with flat array-like objects works", {
    set.seed(107)

    # Array-like objects (not SE)
    mat1 <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
    mat2 <- matrix(rpois(60 * 20, 10), nrow = 60, ncol = 20)

    rownames(mat1) <- paste0("Feature1_", 1:50)
    rownames(mat2) <- paste0("Feature2_", 1:60)
    colnames(mat1) <- colnames(mat2) <- paste0("Sample", 1:20)

    mae <- MultiAssayExperiment(
        experiments = ExperimentList(Dataset1 = mat1, Dataset2 = mat2)
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(mae, tmpdir)
    mae2 <- readParquet(tmpdir)

    # Check structure
    expect_s4_class(mae2, "MultiAssayExperiment")
    expect_identical(names(experiments(mae2)), c("Dataset1", "Dataset2"))
    expect_identical(length(experiments(mae2)), 2L)

    # Array-like objects come back as DuckDBMatrix
    expect_s4_class(mae2[["Dataset1"]], "DuckDBMatrix")
    expect_s4_class(mae2[["Dataset2"]], "DuckDBMatrix")

    # Check matrices with proper dimnames names
    mat1_expected <- mat1
    names(dimnames(mat1_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(mae2[["Dataset1"]], mat1_expected)

    mat2_expected <- mat2
    names(dimnames(mat2_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(mae2[["Dataset2"]], mat2_expected)

    unlink(tmpdir, recursive = TRUE)
})

# ==============================================================================
# Unbound Metadata Tests
# ==============================================================================

test_that("Unbound metadata with serializable S4 objects works", {
    set.seed(108)
    ncells <- 20L
    ngenes <- 30L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    qc_df <- DataFrame(
        metric = c("n_genes", "n_counts"),
        threshold = c(200, 1000)
    )

    # Create SE with complex metadata
    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = DataFrame(cluster = sample(1:2, ncells, replace = TRUE)),
        metadata = list(
            version = "1.0",
            qc_stats = qc_df
        )
    )

    # Round-trip
    tmpdir <- tempfile()
    writeParquet(se, tmpdir)
    se2 <- readParquet(tmpdir)

    # Check structure
    expect_identical(dim(se2), dim(se))
    counts_expected <- assay(se, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(se2, "counts"), counts_expected)
    checkDuckDBDataFrame(colData(se2), as.data.frame(colData(se)))

    # Check metadata
    expect_identical(metadata(se2)$version, "1.0")
    expect_s4_class(metadata(se2)$qc_stats, "DuckDBDataFrame")
    checkDuckDBDataFrame(metadata(se2)$qc_stats, as.data.frame(qc_df))

    unlink(tmpdir, recursive = TRUE)
})

test_that("Nested metadata serializes JSON leaves and tabular sidecars", {
    set.seed(110)
    ncells <- 20L
    ngenes <- 30L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = DataFrame(cluster = sample(1:2, ncells, replace = TRUE)),
        metadata = list(
            bulk_labels_colors = array(rep("#FF0000", 10L)),
            louvain = list(params = list(resolution = 1.0)),
            pca = list(variance = rnorm(50L), variance_ratio = rnorm(50L)),
            rank_genes_groups = list(
                names = matrix(letters[1:20], 10L, 2L),
                scores = matrix(rnorm(20L), 10L, 2L),
                params = list(groupby = "bulk_labels")
            ),
            qc_stats = DataFrame(metric = "n_genes", threshold = 200L)
        )
    )

    tmpdir <- tempfile()
    writeParquet(se, tmpdir)
    pkg <- jsonlite::fromJSON(file.path(tmpdir, "datapackage.json"),
                              simplifyVector = FALSE)

    expect_type(pkg$annotations$bulk_labels_colors, "list")
    expect_null(pkg$annotations$bulk_labels_colors$`__type__`)
    expect_equal(pkg$annotations$louvain$params$resolution, 1.0)
    expect_length(pkg$annotations$pca$variance, 50L)
    expect_identical(pkg$annotations$rank_genes_groups$`__type__`,
                     "nested_mapping")
    expect_identical(pkg$annotations$rank_genes_groups$names$`__type__`,
                     "parquet_ref")
    expect_identical(pkg$annotations$rank_genes_groups$params$groupby,
                     "bulk_labels")
    expect_identical(pkg$annotations$qc_stats$`__type__`, "parquet_ref")

    unbound_names <- vapply(
        Filter(function(r) identical(r$dimension, "unbound"), pkg$resources),
        function(r) r$name,
        character(1L)
    )
    expect_true("rank_genes_groups__names" %in% unbound_names)
    expect_true("rank_genes_groups__scores" %in% unbound_names)
    expect_true("qc_stats" %in% unbound_names)
    expect_false("bulk_labels_colors" %in% unbound_names)
    expect_true(file.exists(file.path(tmpdir, "unbound_rank_genes_groups__names",
                                      "part-0.parquet")))

    se2 <- readParquet(tmpdir)
    expect_equal(metadata(se2)$bulk_labels_colors,
                 as.vector(metadata(se)$bulk_labels_colors))
    expect_equal(metadata(se2)$louvain, metadata(se)$louvain)
    expect_equal(as.numeric(metadata(se2)$pca$variance),
                 as.numeric(metadata(se)$pca$variance), tolerance = 1e-4)
    expect_identical(metadata(se2)$rank_genes_groups$params,
                     metadata(se)$rank_genes_groups$params)
    expect_s4_class(metadata(se2)$rank_genes_groups$names, "DuckDBDataFrame")
    checkDuckDBDataFrame(metadata(se2)$rank_genes_groups$names,
                         as.data.frame(metadata(se)$rank_genes_groups$names))
    checkDuckDBDataFrame(metadata(se2)$rank_genes_groups$scores,
                         as.data.frame(metadata(se)$rank_genes_groups$scores))
    expect_s4_class(metadata(se2)$qc_stats, "DuckDBDataFrame")
    checkDuckDBDataFrame(metadata(se2)$qc_stats,
                         as.data.frame(metadata(se)$qc_stats))

    unlink(tmpdir, recursive = TRUE)
})

test_that("Metadata with package_version and POSIXt serializes to JSON", {
    set.seed(111)
    se <- SummarizedExperiment(
        assays = list(counts = matrix(1:12, 3, 4)),
        metadata = list(
            package_version = packageVersion("BiocDuckDB"),
            creation_date = as.POSIXct("2026-01-01 12:00:00", tz = "UTC")
        )
    )
    tmpdir <- tempfile()
    expect_silent(writeParquet(se, tmpdir))
    pkg <- jsonlite::fromJSON(file.path(tmpdir, "datapackage.json"))
    expect_type(pkg$annotations$package_version, "character")
    expect_type(pkg$annotations$creation_date, "character")
    se2 <- readParquet(tmpdir)
    expect_equal(metadata(se2)$package_version, as.character(packageVersion("BiocDuckDB")))
    unlink(tmpdir, recursive = TRUE)
})

test_that("Unbound metadata with non-serializable objects is skipped with warning", {
    # Load airway which has MIAME metadata
    data(airway, package = "airway")

    # Check that MIAME metadata exists
    expect_s4_class(metadata(airway)[[1]], "MIAME")

    # Round-trip should succeed (not error) but MIAME is skipped
    tmpdir <- tempfile()
    expect_warning(writeParquet(airway, tmpdir), "Skipping unsupported metadata")

    airway2 <- readParquet(tmpdir)

    # MIAME should not be in reconstructed object (silently dropped)
    expect_false("MIAME" %in% vapply(metadata(airway2), function(x) class(x)[1], character(1L)))

    # But the rest of the object should be intact
    expect_identical(dim(airway2), dim(airway))
    airway_counts_expected <- assay(airway, "counts")
    names(dimnames(airway_counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(airway2, "counts"), airway_counts_expected)

    # colData factors now round-trip, so compare factor-vs-factor
    coldata_expected <- as.data.frame(colData(airway), stringsAsFactors = FALSE)
    checkDuckDBDataFrame(colData(airway2), coldata_expected)

    # Check rowRanges (basic structure only - GRangesList metadata has known limitations)
    expect_s4_class(rowRanges(airway2), "DuckDBGRangesList")
    expect_identical(length(rowRanges(airway2)), length(rowRanges(airway)))
    expect_identical(names(rowRanges(airway2)), names(rowRanges(airway)))

    unlink(tmpdir, recursive = TRUE)
})

test_that("factor and ordered factor colData columns survive round-trip", {
    ncells <- 12L
    counts <- matrix(rpois(20L * ncells, 5), nrow = 20L)
    rownames(counts) <- paste0("Gene", seq_len(20L))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        colData = DataFrame(
            group = factor(rep(c("ctrl", "treat", "ctrl"), length.out = ncells),
                           levels = c("ctrl", "treat")),
            dose = factor(rep(c("low", "high", "mid"), length.out = ncells),
                          levels = c("low", "mid", "high"), ordered = TRUE),
            label = rep(c("x", "y"), length.out = ncells)  # plain character
        )
    )

    tmpdir <- tempfile()
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
    writeParquet(se, tmpdir)
    se2 <- readParquet(tmpdir)

    cd <- as.data.frame(colData(se2))
    expect_true(is.factor(cd$group))
    expect_false(is.ordered(cd$group))
    expect_identical(levels(cd$group), c("ctrl", "treat"))

    expect_true(is.ordered(cd$dose))
    expect_identical(levels(cd$dose), c("low", "mid", "high"))

    expect_type(cd$label, "character")

    # single-column materialization path (as.vector,DuckDBColumn)
    expect_true(is.ordered(as.vector(colData(se2)[["dose"]])))
})

# ==============================================================================
# Model Dispatch Tests
# ==============================================================================

test_that("NULL model returns SimpleList", {
    set.seed(109)

    # Create a minimal datapackage.json with no model field
    tmpdir <- tempfile()
    dir.create(tmpdir)

    # Write a simple DataFrame resource
    df <- DataFrame(x = 1:10, y = letters[1:10])
    df_path <- file.path(tmpdir, "data")
    writeParquet(df, df_path, dimension = "unbound")

    # Manually create datapackage.json without model field
    pkg <- list(
        resources = list(
            list(
                name = "data",
                path = "data",
                dimension = "unbound",
                layout = "data_frame",
                format = "parquet",
                mediatype = "application/vnd.apache.parquet"
            )
        )
    )
    jsonlite::write_json(pkg, file.path(tmpdir, "datapackage.json"), 
                        auto_unbox = TRUE)

    # Read without model
    result <- readParquet(tmpdir)

    # Should return SimpleList
    expect_s4_class(result, "SimpleList")
    expect_identical(names(result), "data")
    expect_s4_class(result[[1]], "DuckDBDataFrame")

    unlink(tmpdir, recursive = TRUE)
})

test_that("Model dispatch routes to correct readers", {
    # SummarizedExperiment
    se <- SummarizedExperiment(assays = list(counts = matrix(1:12, 3, 4)))
    tmpdir1 <- tempfile()
    writeParquet(se, tmpdir1)
    pkg1 <- jsonlite::read_json(file.path(tmpdir1, "datapackage.json"))
    expect_identical(pkg1$model, "summarized_experiment")
    se2 <- readParquet(tmpdir1)
    expect_s4_class(se2, "SummarizedExperiment")
    unlink(tmpdir1, recursive = TRUE)

    # SingleCellExperiment
    sce <- SingleCellExperiment(assays = list(counts = matrix(1:12, 3, 4)))
    tmpdir2 <- tempfile()
    writeParquet(sce, tmpdir2)
    pkg2 <- jsonlite::read_json(file.path(tmpdir2, "datapackage.json"))
    expect_identical(pkg2$model, "single_cell_experiment")
    sce2 <- readParquet(tmpdir2)
    expect_s4_class(sce2, "SingleCellExperiment")
    unlink(tmpdir2, recursive = TRUE)

    # Note: ExperimentList tests skipped - List serialization needs explicit layout
})

# ==============================================================================
# Dimension and Layout Tests
# ==============================================================================

test_that("Dimension and layout fields are correctly written", {
    set.seed(110)
    ncells <- 15L
    ngenes <- 20L

    counts <- matrix(rpois(ngenes * ncells, 5), nrow = ngenes, ncol = ncells)
    rownames(counts) <- paste0("Gene", seq_len(ngenes))
    colnames(counts) <- paste0("Cell", seq_len(ncells))

    sce <- SingleCellExperiment(
        assays = list(counts = counts),
        rowData = DataFrame(gene_type = sample(c("A", "B"), ngenes, replace = TRUE)),
        colData = DataFrame(cluster = sample(1:2, ncells, replace = TRUE))
    )

    tmpdir <- tempfile()
    writeParquet(sce, tmpdir)

    # Read datapackage.json
    pkg <- jsonlite::read_json(file.path(tmpdir, "datapackage.json"),
                              simplifyVector = TRUE)

    # Check model
    expect_identical(pkg$model, "single_cell_experiment")

    # Check resource dimensions and layouts
    resources <- pkg$resources
    if (nrow(resources) > 0) {
        # resources is a data.frame when simplifyVector = TRUE
        rownames(resources) <- resources$name

        # Check that dimension and layout fields exist and are non-empty
        expect_true("features" %in% resources$name)
        expect_true(!is.na(resources["features", "dimension"]))
        expect_true(!is.na(resources["features", "layout"]))
        expect_identical(resources["features", "dimension"], "feature")

        expect_true("samples" %in% resources$name)
        expect_identical(resources["samples", "dimension"], "sample")
        expect_identical(resources["samples", "layout"], "data_frame")

        expect_true("counts" %in% resources$name)
        expect_identical(resources["counts", "dimension"], "crossed")
        expect_identical(resources["counts", "layout"], "coord_array")
    }

    unlink(tmpdir, recursive = TRUE)
})

cat("All readParquet/writeParquet tests defined.\n")
