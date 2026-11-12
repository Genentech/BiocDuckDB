test_that("writeParquet falls back to snappy when zstd is unavailable", {
    skip_if_not_installed("arrow")
    skip_if_not_installed("nanoparquet")
    testthat::local_mocked_bindings(codec_is_available = function(type) FALSE, .package = "arrow")

    path <- tempfile()
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)

    df <- data.frame(value = 1:5)
    writeParquet(df, path, indexcol = "__sample__", keycol = "__name__",
                part = 0L, name = "samples", dimension = "sample")

    pq <- file.path(path, "part-0.parquet")
    got <- arrow::read_parquet(pq)
    expect_equal(got$value, 1:5)
    codecs <- unique(nanoparquet::read_parquet_metadata(pq)$column_chunks$codec)
    expect_identical(codecs, "SNAPPY")
})
