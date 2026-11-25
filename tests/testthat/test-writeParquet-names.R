# writeParquet() does not rename a column behind the back of the caller. DuckDB,
# which reads these files, binds a column name without regard to case, so a file
# with both Sex and sex gives the values of Sex for either; and a column named
# like the index or the key used to be renamed with a suffix.

.parquet_columns <- function(dir) {
    pq <- list.files(dir, pattern = "[.]parquet$", recursive = TRUE, full.names = TRUE)
    names(arrow::read_parquet(pq[[1L]], as_data_frame = FALSE)$schema)
}

test_that("writeParquet refuses columns whose names are the same up to case", {
    mem <- S4Vectors::DataFrame(Sex = c("M", "F", "M"),
                                sex = c("male", "female", "male"),
                                row.names = c("a", "b", "c"), check.names = FALSE)
    expect_identical(colnames(mem), c("Sex", "sex"))
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE))
    expect_error(writeParquet(mem, out), "'Sex', 'sex'")
    expect_error(writeParquet(as.data.frame(mem), tempfile()), "same name, ignoring case")
})

test_that("writeParquet refuses repeated column names instead of renaming them", {
    dup <- S4Vectors::DataFrame(x = c(1, 2), x = c(10, 20), check.names = FALSE)
    expect_identical(colnames(dup), c("x", "x"))
    expect_error(writeParquet(dup, tempfile()), "same name, ignoring case")
})

test_that("writeParquet refuses a column named like the index or the key", {
    # it used to be renamed __name___1, silently
    mem <- S4Vectors::DataFrame(`__name__` = c("U", "V"), v = 1:2,
                                row.names = c("a", "b"), check.names = FALSE)
    expect_error(writeParquet(mem, tempfile()), "'__name__', '__name__'")
    idx <- S4Vectors::DataFrame(`__index__` = c(5L, 6L), v = 1:2, check.names = FALSE)
    expect_error(writeParquet(idx, tempfile()), "'__index__', '__index__'")
})

test_that("writeParquet keeps names that are not valid R names, and names that differ by more than case", {
    mem <- S4Vectors::DataFrame(`CD3-positive` = c(1, 2), `n cells` = c(10, 20),
                                `2nd` = c("a", "b"), Sex = c("M", "F"), Sex2 = c("m", "f"),
                                row.names = c("r1", "r2"), check.names = FALSE)
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE))
    writeParquet(mem, out)
    expect_identical(.parquet_columns(out),
                     c("__index__", "__name__", "CD3-positive", "n cells", "2nd", "Sex", "Sex2"))
})

test_that("a lazily read table with names that clash up to case rewrites the column that it exposes", {
    d <- data.frame(id = c("a", "b", "c"), Sex = c("M", "F", "M"),
                    sex = c("male", "female", "male"), check.names = FALSE)
    src <- tempfile(fileext = ".parquet")
    arrow::write_parquet(d, src)
    on.exit(unlink(src))
    both <- DuckDBDataFrame::DuckDBDataFrame(src, datacols = c("Sex", "sex"), keycol = "id")
    expect_error(writeParquet(both, tempfile()), "'Sex', 'sex'")

    # one of the two columns is written with its own values
    out <- tempfile()
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    only <- DuckDBDataFrame::DuckDBDataFrame(src, datacols = "sex", keycol = "id")
    writeParquet(only, out)
    pq <- list.files(out, pattern = "[.]parquet$", recursive = TRUE, full.names = TRUE)[[1L]]
    written <- as.data.frame(arrow::read_parquet(pq))
    expect_identical(written$sex, d$sex)
})
