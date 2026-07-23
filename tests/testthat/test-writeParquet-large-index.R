# Regression: append offsets and __index__ values above the 32-bit ceiling (a
# resource with more than ~2.1e9 rows). `index_max` declares the index range so
# the __index__ column is typed wide enough (int64) up front -- part 0 included
# -- and every streamed part shares one type, instead of narrowing part 0 to
# int32 and overflowing later parts on cast. See NEWS 0.99.3.
# Run: library(BiocDuckDB); library(testthat); source("test-writeParquet-large-index.R")

library(arrow)

.indexArrowType <- function(f) {
    ParquetFileReader$create(f)$GetSchema()$GetFieldByName("__index__")$type$ToString()
}

test_that("index_max streams a > 2^31 index as a consistent int64 column", {
    dir <- tempfile()
    # Part 0 at offset 0, then a part whose offset crosses the 32-bit ceiling.
    writeParquet(data.frame(v = 1:5), dir, indexcol = "__index__", keycol = NULL,
                 dimension = "sample", layout = "data_frame",
                 offset = 0, part = 0L, part_digits = 2L, append = FALSE,
                 index_max = Inf)
    writeParquet(data.frame(v = 11:15), dir, indexcol = "__index__", keycol = NULL,
                 dimension = "sample", layout = "data_frame",
                 offset = 3e9, part = 1L, part_digits = 2L, append = TRUE,
                 index_max = Inf)

    files <- sort(list.files(dir, pattern = "parquet$", recursive = TRUE,
                             full.names = TRUE))
    expect_length(files, 2L)
    # Both parts, and the unified dataset, are int64 (schema-consistent).
    expect_true(all(vapply(files, .indexArrowType, character(1L)) == "int64"))
    ds <- open_dataset(dir)
    expect_identical(ds$schema$GetFieldByName("__index__")$type$ToString(), "int64")

    # The index values above 2^31 are stored exactly, with no overflow to NA.
    idx <- sort(as.data.frame(ds)[["__index__"]])
    expect_false(anyNA(idx))
    expect_equal(idx[6:10],
                 c(3000000001, 3000000002, 3000000003, 3000000004, 3000000005))
    unlink(dir, recursive = TRUE)
})

test_that("without index_max a small index still narrows (no regression)", {
    dir <- tempfile()
    writeParquet(data.frame(v = 1:5), dir, indexcol = "__index__", keycol = NULL,
                 dimension = "sample", layout = "data_frame")
    f <- list.files(dir, pattern = "parquet$", recursive = TRUE,
                    full.names = TRUE)[1L]
    expect_identical(.indexArrowType(f), "uint8")
    unlink(dir, recursive = TRUE)
})

test_that("a > 2^31 index without index_max fails loudly (not silent float64)", {
    dir <- tempfile()
    expect_error(
        writeParquet(data.frame(v = 1:5), dir, indexcol = "__index__",
                     keycol = NULL, dimension = "sample", layout = "data_frame",
                     offset = 3e9, part = 0L, part_digits = 2L, append = FALSE),
        "index_max")
    unlink(dir, recursive = TRUE)
})

test_that("cross-path append is schema-consistent", {
    skip_if_not_installed("arrow")
    dir <- tempfile()
    # Part 0 via the in-memory data.frame path, typed int64 up front.
    writeParquet(data.frame(v = 1:5), dir, indexcol = "__index__", keycol = NULL,
                 dimension = "sample", layout = "data_frame",
                 offset = 0, part = 0L, part_digits = 2L, append = FALSE,
                 index_max = Inf)
    # Part 1 via the lazy DuckDBTable path; must pin __index__ to part 0's int64
    # (not default to BIGINT independently) so the resource stays readable.
    src <- tempfile(fileext = ".parquet"); on.exit(unlink(src), add = TRUE)
    arrow::write_parquet(data.frame(v = 11:15), src)
    ddf <- DuckDBDataFrame::DuckDBDataFrame(src)
    writeParquet(ddf, dir, indexcol = "__index__", keycol = NULL,
                 dimension = "sample", layout = "data_frame",
                 offset = 5, part = 1L, part_digits = 2L, append = TRUE)

    files <- sort(list.files(dir, pattern = "parquet$", recursive = TRUE,
                             full.names = TRUE))
    expect_length(files, 2L)
    types <- vapply(files, function(p)
        arrow::ParquetFileReader$create(p)$GetSchema()$
            GetFieldByName("__index__")$type$ToString(), character(1L))
    expect_true(all(types == "int64"))
    unlink(dir, recursive = TRUE)
})
