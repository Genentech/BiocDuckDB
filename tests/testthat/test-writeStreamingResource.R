# Tests for writeStreamingResource(): the public streaming/append writer that
# owns the offset / part / append / part_digits / dimtbl-slice / index_max
# bookkeeping over writeParquet (review finding F18).

library(S4Vectors)

test_that("streams multiple blocks into one contiguous multi-part resource", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    chunks <- list(
        data.frame(v = 1:3, g = c("x", "y", "z"), stringsAsFactors = FALSE),
        data.frame(v = 4:5, g = c("p", "q"), stringsAsFactors = FALSE))
    res <- writeStreamingResource(
        function(i) if (i <= length(chunks)) chunks[[i]] else NULL,
        path = file.path(tf, "tbl"), dimension = "unbound", keycol = NULL,
        index_max = 5, expected_rows = 5, name = "tbl")

    expect_type(res, "list")
    expect_identical(res[[1L]][["n_rows"]], 5)
    expect_length(list.files(file.path(tf, "tbl"), pattern = "\\.parquet$"), 2L)

    # Parts read back as one table with a contiguous 1..5 index and intact values.
    tbl <- DuckDBDataFrame::DuckDBDataFrame(file.path(tf, "tbl"),
                                            keycol = "__index__")
    df <- as.data.frame(tbl)
    expect_identical(sort(as.integer(rownames(df))), 1:5)
    expect_identical(sort(df$v), 1:5)
    expect_identical(sort(df$g), c("p", "q", "x", "y", "z"))
})

test_that("empty and leading-empty streams are handled", {
    expect_null(writeStreamingResource(
        function(i) NULL, path = file.path(tempfile(), "t"),
        dimension = "unbound", keycol = NULL))

    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    # a leading zero-row block must not consume part 0
    res <- writeStreamingResource(
        function(i) switch(i,
            data.frame(v = integer(0)),   # empty leading block
            data.frame(v = 1:4),          # -> becomes part 0
            NULL),
        path = file.path(tf, "t"), dimension = "unbound", keycol = NULL,
        index_max = 4)
    expect_identical(res[[1L]][["n_rows"]], 4)
    expect_identical(list.files(file.path(tf, "t"), pattern = "\\.parquet$"),
                     "part-000000.parquet")
})

test_that("warns on the narrowing floor: small part 0, more parts, no index_max", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    expect_warning(
        writeStreamingResource(
            function(i) if (i <= 2) data.frame(v = 1:10) else NULL,
            path = file.path(tf, "t"), dimension = "unbound", keycol = NULL),
        "index_max")

    # ... but not when index_max is supplied
    tf2 <- tempfile()
    on.exit(unlink(tf2, recursive = TRUE), add = TRUE)
    expect_no_warning(
        writeStreamingResource(
            function(i) if (i <= 2) data.frame(v = 1:10) else NULL,
            path = file.path(tf2, "t"), dimension = "unbound", keycol = NULL,
            index_max = 20))
})

test_that("warns when the streamed row count misses expected_rows", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    expect_warning(
        writeStreamingResource(
            function(i) if (i == 1L) data.frame(v = 1:2) else NULL,
            path = file.path(tf, "t"), dimension = "unbound", keycol = NULL,
            expected_rows = 5),
        "expected")
})

test_that("fails fatally when the stream and dimtbl row counts disagree", {
    dt <- DataFrame(grp = c("a", "a", "b"), row.names = as.character(1:3))
    expect_error(
        writeStreamingResource(
            function(i) if (i == 1L) data.frame(v = 1:2) else NULL,
            path = file.path(tempfile(), "t"), dimension = "sample",
            keycol = NULL, dimtbl = dt, index_max = 3),
        "misaligned")
})

test_that("a non-data.frame block is coerced", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    res <- writeStreamingResource(
        function(i) if (i == 1L) DataFrame(v = 1:3) else NULL,
        path = file.path(tf, "t"), dimension = "unbound", keycol = NULL,
        index_max = 3)
    expect_identical(res[[1L]][["n_rows"]], 3)
})
