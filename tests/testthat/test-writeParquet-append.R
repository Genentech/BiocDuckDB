test_that("writeParquet flat parts use offset for index column", {
    path <- tempfile()
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)

    df1 <- data.frame(value = 1:3)
    rownames(df1) <- c("a", "b", "c")
    df2 <- data.frame(value = 4:5)
    rownames(df2) <- c("d", "e")

    pkg <- list(resources = list())
    res1 <- writeParquet(df1, path, indexcol = "__sample__", keycol = "__name__",
                         part = 0L, name = "samples", dimension = "sample")
    expect_length(res1, 1L)
    pkg$resources <- c(pkg$resources, res1)

    res2 <- writeParquet(df2, path, indexcol = "__sample__", keycol = "__name__",
                         offset = 3L, part = 1L, append = TRUE,
                         name = "samples", dimension = "sample")
    expect_null(res2)
    expect_length(pkg$resources, 1L)

  pq1 <- arrow::read_parquet(file.path(path, "part-0.parquet"))
  pq2 <- arrow::read_parquet(file.path(path, "part-1.parquet"))
  expect_equal(pq1$`__sample__`, 1:3)
  expect_equal(pq2$`__sample__`, 4:5)
})

test_that("writeParquet part_digits zero-pads filenames", {
    path <- tempfile()
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)

    df <- data.frame(x = 1L)
    writeParquet(df, path, indexcol = NULL, part = 0L, part_digits = 2L)
    expect_true(file.exists(file.path(path, "part-00.parquet")))
})

test_that("writeParquet append refuses duplicate part", {
    path <- tempfile()
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)

    df <- data.frame(x = 1L)
    writeParquet(df, path, indexcol = NULL, part = 0L)
    expect_error(
        writeParquet(df, path, indexcol = NULL, part = 0L, append = TRUE),
        "already exists"
    )
})

test_that("writeParquet append requires part", {
    path <- tempfile()
    dir.create(path)
    on.exit(unlink(path, recursive = TRUE), add = TRUE)

    df <- data.frame(x = 1L)
    expect_error(
        writeParquet(df, path, append = TRUE),
        "requires 'part'"
    )
})
