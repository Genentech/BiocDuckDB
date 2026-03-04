# Tests that spatialOverlaps and spatialMatch dispatch correctly for
# DuckDBDataFrame and produce results identical to the DataFrame (sf) path.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-spatial-dispatch.R")

skip_if_not_installed("MultiAssaySpatialExperiment")
skip_if_not_installed("sf")

library(sf)
library(MultiAssaySpatialExperiment)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Test data
###

.make_points_ddb <- function() {
    df <- data.frame(
        id = 1:6,
        x = c(1, 5, 10, 50, 100, 3),
        y = c(1, 5, 10, 50, 100, 3)
    )
    path <- tempfile(fileext = ".parquet")
    arrow::write_parquet(df, path)
    DuckDBDataFrame(path, datacols = c("x", "y"), keycol = "id")
}

.make_points_df <- function() {
    DataFrame(x = c(1, 5, 10, 50, 100, 3),
              y = c(1, 5, 10, 50, 100, 3),
              row.names = as.character(1:6))
}

.make_shapes_df <- function() {
    wkt <- c("POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))",
             "POLYGON((7 7, 12 7, 12 12, 7 12, 7 7))",
             "POLYGON((40 40, 60 40, 60 60, 40 60, 40 40))")
    DataFrame(geometry = st_as_sfc(wkt))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialOverlaps
###

test_that("spatialOverlaps dispatches to DuckDBDataFrame method (coords mode)", {
    ddb <- .make_points_ddb()
    mem <- .make_points_df()
    polygon <- st_as_sfc("POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))")

    ddb_result <- spatialOverlaps(ddb, polygon, coords = c("x", "y"))
    mem_result <- spatialOverlaps(mem, polygon, coords = c("x", "y"))

    expect_true(is.logical(ddb_result))
    expect_identical(length(ddb_result), nrow(ddb))
    expect_identical(unname(ddb_result), unname(mem_result))
})

test_that("spatialOverlaps DuckDB result matches sf for various polygons", {
    ddb <- .make_points_ddb()
    mem <- .make_points_df()

    big_poly <- st_as_sfc("POLYGON((0 0, 200 0, 200 200, 0 200, 0 0))")
    expect_identical(
        unname(spatialOverlaps(ddb, big_poly, coords = c("x", "y"))),
        unname(spatialOverlaps(mem, big_poly, coords = c("x", "y"))))

    small_poly <- st_as_sfc("POLYGON((0 0, 2 0, 2 2, 0 2, 0 0))")
    expect_identical(
        unname(spatialOverlaps(ddb, small_poly, coords = c("x", "y"))),
        unname(spatialOverlaps(mem, small_poly, coords = c("x", "y"))))

    no_hit <- st_as_sfc("POLYGON((500 500, 600 500, 600 600, 500 600, 500 500))")
    expect_identical(
        unname(spatialOverlaps(ddb, no_hit, coords = c("x", "y"))),
        unname(spatialOverlaps(mem, no_hit, coords = c("x", "y"))))
})

test_that("spatialOverlaps accepts WKT string for y", {
    ddb <- .make_points_ddb()
    polygon_wkt <- "POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))"
    result <- spatialOverlaps(ddb, polygon_wkt, coords = c("x", "y"))
    expect_true(is.logical(result))
    expect_true(any(result))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialMatch
###

test_that("spatialMatch dispatches to DuckDBDataFrame method", {
    ddb <- .make_points_ddb()
    mem <- .make_points_df()
    shapes <- .make_shapes_df()

    ddb_result <- spatialMatch(ddb, shapes, coords = c("x", "y"))
    mem_result <- spatialMatch(mem, shapes, coords = c("x", "y"))

    expect_true(is.integer(ddb_result))
    expect_identical(length(ddb_result), nrow(ddb))
    expect_identical(ddb_result, mem_result)
})

test_that("spatialMatch returns NA for unmatched points", {
    ddb <- .make_points_ddb()
    shapes <- .make_shapes_df()
    result <- spatialMatch(ddb, shapes, coords = c("x", "y"))
    expect_true(is.na(result[5L]))
})

test_that("spatialMatch assigns correct region indices", {
    ddb <- .make_points_ddb()
    shapes <- .make_shapes_df()
    result <- spatialMatch(ddb, shapes, coords = c("x", "y"))
    expect_identical(result[1L], 1L)
    expect_identical(result[3L], 2L)
    expect_identical(result[4L], 3L)
})
