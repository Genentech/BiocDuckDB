# Tests that spatialOverlaps, spatialMatch, and spatialJoin dispatch correctly
# for DuckDBDataFrame and produce results identical to the DataFrame (sf) path.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-spatial-dispatch.R")

skip_if_not_installed("MultiAssaySpatialExperiment")
skip_if_not_installed("sf")
skip_if_not_installed("DuckDBSpatial")

# The spatial predicates need DuckDB's spatial extension to be loadable; skip
# where it cannot be obtained (offline / firewalled extension repository).
skip_if_not(
    tryCatch({
        DBI::dbGetQuery(DuckDBDataFrame::acquireDuckDBConn(),
            "SELECT ST_Area(ST_GeomFromText('POLYGON((0 0,0 1,1 1,1 0,0 0))'))")
        TRUE
    }, error = function(e) FALSE),
    "DuckDB spatial extension not available")

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

.make_shapes_sf <- function() {
    wkt <- c("POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))",
             "POLYGON((7 7, 12 7, 12 12, 7 12, 7 7))",
             "POLYGON((40 40, 60 40, 60 60, 40 60, 40 40))")
    st_sf(id = seq_along(wkt), geometry = st_as_sfc(wkt))
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

# 'sf' is a base S3 class (c("sf", "data.frame")), not an S4Vectors DataFrame,
# so the DataFrame-signature methods above never apply to it; a dedicated
# 'sf'-signature method is required for S4 dispatch to find it at all.
test_that("spatialMatch dispatches to the DuckDB-backed method for a plain sf table", {
    ddb <- .make_points_ddb()
    mem_result <- spatialMatch(.make_points_df(), .make_shapes_df(), coords = c("x", "y"))
    sf_result <- spatialMatch(ddb, .make_shapes_sf(), coords = c("x", "y"))
    expect_identical(sf_result, mem_result)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialJoin
###

# A DuckDBDataFrame joined against a plain in-memory S4Vectors DataFrame (on
# either side) must dispatch to the DuckDBSpatial-backed method, not silently
# fall through to MultiAssaySpatialExperiment's own in-memory
# DataFrame,DataFrame method (which forces both sides fully into memory).
test_that("spatialJoin dispatches to the DuckDB-backed method for a mixed DuckDBDataFrame/DataFrame pair", {
    spatial_path <- system.file("extdata", "spatial", package = "DuckDBSpatial")
    shapes_ddb <- DuckDBDataFrame(spatial_path)
    shapes_ddb <- shapes_ddb[!is.na(shapes_ddb$type), ]
    shapes_mem <- .make_shapes_df()

    forward <- spatialJoin(shapes_ddb, shapes_mem)
    expect_s4_class(forward, "DuckDBDataFrame")
    expect_true(nrow(forward) > 0L)

    backward <- spatialJoin(shapes_mem, shapes_ddb)
    expect_s4_class(backward, "DuckDBDataFrame")
    expect_identical(nrow(forward), nrow(backward))
})

# Same gap, but for 'sf' specifically: it's a base S3 class, not an
# S4Vectors DataFrame, so it needs its own signature for S4 dispatch to
# find a method at all (it previously errored with "unable to find an
# inherited method" rather than falling through to an in-memory default).
test_that("spatialJoin dispatches to the DuckDB-backed method for a mixed DuckDBDataFrame/sf pair", {
    spatial_path <- system.file("extdata", "spatial", package = "DuckDBSpatial")
    shapes_ddb <- DuckDBDataFrame(spatial_path)
    shapes_ddb <- shapes_ddb[!is.na(shapes_ddb$type), ]
    shapes_sf <- .make_shapes_sf()

    forward <- spatialJoin(shapes_ddb, shapes_sf)
    expect_s4_class(forward, "DuckDBDataFrame")
    expect_true(nrow(forward) > 0L)

    backward <- spatialJoin(shapes_sf, shapes_ddb)
    expect_s4_class(backward, "DuckDBDataFrame")
    expect_identical(nrow(forward), nrow(backward))
})
