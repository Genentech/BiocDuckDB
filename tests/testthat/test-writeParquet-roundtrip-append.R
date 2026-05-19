# Round-trip tests for incremental writeParquet (flat multi-part + hive append)

library(SummarizedExperiment)

.equalTicks <- function(n, k) {
    as.integer(round(seq(0, n, length.out = k + 1L)))[-1L]
}

.groupFromTicks <- function(indices, ticks) {
    findInterval(indices - 1L, ticks) + 1L
}

.featureDimtbl <- function(n_features, feature_ticks, indexcol = "__feature__",
                           grid_suffix = "group__") {
    df <- data.frame(x = as.integer(.groupFromTicks(seq_len(n_features),
                                                    feature_ticks)))
    names(df) <- paste0(indexcol, grid_suffix)
    df
}

.sampleDimtbl <- function(n_samples, sample_ticks, group_offset,
                          indexcol = "__sample__", grid_suffix = "group__") {
    local <- .groupFromTicks(seq_len(n_samples), sample_ticks)
    df <- data.frame(x = as.integer(local + group_offset))
    names(df) <- paste0(indexcol, grid_suffix)
    df
}

.assayCountsResource <- function(path, indexcols, datacol, name = "counts") {
    fields <- c(lapply(indexcols, function(nm) list(name = nm)),
                list(list(name = datacol)))
    foreignKeys <- list(
        list(fields = indexcols[1L],
             reference = list(fields = indexcols[1L], resource = "features")),
        list(fields = indexcols[2L],
             reference = list(fields = indexcols[2L], resource = "samples"))
    )
    list(name = name,
         path = basename(path),
         dimension = "crossed",
         layout = "coord_array",
         format = "parquet",
         mediatype = "application/vnd.apache.parquet",
         schema = list(fields = fields, foreignKeys = foreignKeys))
}

test_that("readParquet reads multi-part flat samples with global index offsets", {
    tmpdir <- tempfile()
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

    indexcol <- "__sample__"
    keycol <- "__name__"
    samples_dir <- file.path(tmpdir, "samples")

    chunks <- list(
        data.frame(batch = "A", x = 1:5, row.names = paste0("S", 1:5)),
        data.frame(batch = "B", x = 6:10, row.names = paste0("S", 6:10)),
        data.frame(batch = "C", x = 11:15, row.names = paste0("S", 11:15))
    )
    offsets <- c(0L, 5L, 10L)

    samples_resource <- NULL
    for (i in seq_along(chunks)) {
        res <- writeParquet(chunks[[i]], samples_dir,
                            indexcol = indexcol, keycol = keycol,
                            offset = offsets[i], part = i - 1L,
                            append = i > 1L,
                            name = "samples", dimension = "sample")
        if (i == 1L) {
            samples_resource <- res[[1L]]
        }
    }

    expect_true(all(file.exists(file.path(samples_dir,
                                          sprintf("part-%d.parquet", 0:2)))))

    pkg <- list(resources = list(samples_resource))
    jsonlite::write_json(pkg, file.path(tmpdir, "datapackage.json"),
                         auto_unbox = TRUE)

    out <- readParquet(tmpdir, model = NULL)
    df <- as.data.frame(out$samples)
    expect_equal(nrow(df), 15L)
    # Global index is the DuckDB key column, not a regular data column
    expect_equal(as.integer(out$samples@keycols[[indexcol]]), seq_len(15L))
    expect_equal(df$batch, rep(c("A", "B", "C"), each = 5L))
    expect_equal(df$x, seq_len(15L))
})

test_that("readParquet round-trips multi-part samples and appended hive counts", {
    set.seed(20260518)
    ngenes <- 6L
    nsamples <- 12L
    slab_size <- 4L
    n_slabs <- nsamples / slab_size

    counts <- matrix(rpois(ngenes * nsamples, 3), nrow = ngenes, ncol = nsamples)
    rownames(counts) <- paste0("G", seq_len(ngenes))
    colnames(counts) <- paste0("S", seq_len(nsamples))

    rowData <- DataFrame(type = rep(c("protein_coding", "lncRNA"), length.out = ngenes))
    rownames(rowData) <- rownames(counts)

    colData <- DataFrame(batch = rep(LETTERS[1:n_slabs], each = slab_size))
    rownames(colData) <- colnames(counts)

    se <- SummarizedExperiment(
        assays = list(counts = counts),
        rowData = rowData,
        colData = colData
    )

    tmpdir <- tempfile()
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

    indexcols <- c("__feature__", "__sample__")
    grid_suffix <- "group__"
    datacol <- "value"

    features_dir <- file.path(tmpdir, "features")
    samples_dir <- file.path(tmpdir, "samples")
    counts_dir <- file.path(tmpdir, "counts")

    feature_ticks <- .equalTicks(ngenes, 2L)
    feature_dimtbl <- .featureDimtbl(ngenes, feature_ticks,
                                     indexcol = indexcols[1L],
                                     grid_suffix = grid_suffix)
    sample_ticks <- .equalTicks(slab_size, 2L)
    max_dim <- c(ngenes, nsamples)

    features_resource <- writeParquet(
        as.data.frame(rowData(se)),
        path = features_dir,
        indexcol = indexcols[1L],
        keycol = "__name__",
        dimtbl = feature_dimtbl,
        part = 0L,
        name = "features",
        dimension = "feature"
    )[[1L]]

    samples_resource <- NULL
    for (i in seq_len(n_slabs)) {
        j <- ((i - 1L) * slab_size + 1L):(i * slab_size)
        sample_dimtbl <- .sampleDimtbl(slab_size, sample_ticks,
                                       group_offset = (i - 1L) * length(sample_ticks),
                                       indexcol = indexcols[2L],
                                       grid_suffix = grid_suffix)
        res <- writeParquet(
            as.data.frame(colData(se)[j, , drop = FALSE]),
            path = samples_dir,
            indexcol = indexcols[2L],
            keycol = "__name__",
            dimtbl = sample_dimtbl,
            offset = (i - 1L) * slab_size,
            part = i - 1L,
            append = i > 1L,
            name = "samples",
            dimension = "sample"
        )
        if (i == 1L) {
            samples_resource <- res[[1L]]
        }
    }

    for (i in seq_len(n_slabs)) {
        j <- ((i - 1L) * slab_size + 1L):(i * slab_size)
        slab <- counts[, j, drop = FALSE]
        grid <- S4Arrays::ArbitraryArrayGrid(list(feature_ticks, sample_ticks))
        if (i == 1L) {
            suppressWarnings(
                writeParquet(slab, counts_dir,
                             indexcols = indexcols,
                             datacol = datacol,
                             grid = grid,
                             grid_suffix = grid_suffix,
                             max_dim = max_dim)
            )
        } else {
            suppressWarnings(
                writeParquet(slab, counts_dir,
                             indexcols = indexcols,
                             datacol = datacol,
                             grid = grid,
                             grid_suffix = grid_suffix,
                             append = TRUE,
                             along = 2L,
                             offset = (i - 1L) * slab_size,
                             group_offset = (i - 1L) * length(sample_ticks))
            )
        }
    }

    expect_true(all(file.exists(file.path(samples_dir,
                                          sprintf("part-%d.parquet", 0:2)))))
    expect_gt(length(list.files(counts_dir, pattern = "\\.parquet$",
                                recursive = TRUE)), 1L)

    counts_resource <- .assayCountsResource(counts_dir, indexcols, datacol)
    package <- list(
        model = "summarized_experiment",
        resources = list(features_resource, samples_resource, counts_resource)
    )
    jsonlite::write_json(package, file.path(tmpdir, "datapackage.json"),
                         auto_unbox = TRUE)

    se2 <- readParquet(tmpdir)

    expect_s4_class(se2, "SummarizedExperiment")
    expect_identical(dim(se2), dim(se))
    checkDuckDBDataFrame(rowData(se2), as.data.frame(rowData(se)))
    checkDuckDBDataFrame(colData(se2), as.data.frame(colData(se)))

    counts_expected <- assay(se, "counts")
    names(dimnames(counts_expected)) <- c("__feature__", "__sample__")
    checkDuckDBMatrix(assay(se2, "counts"), counts_expected)
})
