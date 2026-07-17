# writeParquet(cluster_by=) integration: the clustering key flows from the public generic
# through both the lazy DuckDBDataFrame SQL-COPY path and the materializing data.frame path.
# library(testthat); library(BiocDuckDB); source("test-writeParquet-cluster.R")

.cluster_step_ratio <- function(got, ref) {
    # mean adjacent-row distance of the written order vs the (unclustered) reference
    s <- mean(sqrt(diff(got$x)^2 + diff(got$y)^2))
    b <- mean(sqrt(diff(ref$x)^2 + diff(ref$y)^2))
    s / b
}

.read_one_parquet <- function(dir) {
    f <- list.files(dir, pattern = "parquet$", recursive = TRUE, full.names = TRUE)[1L]
    as.data.frame(arrow::read_parquet(f))
}

test_that("writeParquet(DuckDBDataFrame, cluster_by = zorder()) clusters via the SQL path", {
    skip_if_not_installed("arrow")
    set.seed(1)
    df <- data.frame(x = runif(2500, 0, 100), y = runif(2500, 0, 100))
    src <- tempfile(fileext = ".parquet"); on.exit(unlink(src), add = TRUE)
    arrow::write_parquet(df, src)
    ddf <- DuckDBDataFrame::DuckDBDataFrame(src)

    out <- tempfile()
    writeParquet(ddf, out, cluster_by = DuckDBDataFrame::zorder(c("x", "y")))
    got <- .read_one_parquet(out)
    expect_equal(nrow(got), nrow(df))
    expect_setequal(got$x, df$x)
    expect_lt(.cluster_step_ratio(got, df), 0.5)   # spatially clustered on disk
})

test_that("writeParquet(data.frame, cluster_by = zorder()) clusters via the host path", {
    skip_if_not_installed("arrow")
    set.seed(1)
    df <- data.frame(x = runif(2500, 0, 100), y = runif(2500, 0, 100))
    out <- tempfile()
    writeParquet(df, out, cluster_by = DuckDBDataFrame::zorder(c("x", "y")))
    got <- .read_one_parquet(out)
    expect_equal(nrow(got), nrow(df))
    expect_setequal(got$x, df$x)
    expect_lt(.cluster_step_ratio(got, df), 0.5)
})

test_that("writeParquet without cluster_by is unchanged (no clustering)", {
    skip_if_not_installed("arrow")
    set.seed(1)
    df <- data.frame(x = runif(1500, 0, 100), y = runif(1500, 0, 100))
    src <- tempfile(fileext = ".parquet"); on.exit(unlink(src), add = TRUE)
    arrow::write_parquet(df, src)
    ddf <- DuckDBDataFrame::DuckDBDataFrame(src)
    out <- tempfile()
    writeParquet(ddf, out)   # default __index__ ordering, not spatial
    got <- .read_one_parquet(out)
    expect_equal(nrow(got), nrow(df))
    # not spatially clustered (ratio ~1); this guards that cluster_by is opt-in
    expect_gt(.cluster_step_ratio(got, df), 0.75)
})
