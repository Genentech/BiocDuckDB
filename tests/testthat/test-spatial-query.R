# Cross-element DuckDB query layer over a MASE (MultiAssaySpatialExperiment-query.R):
# on-the-fly views, the spatialMap junction linker, referential-integrity checks,
# coordinate systems, and ST_* cross-element joins with CT-graph alignment.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-spatial-query.R")

# A points-only lazy (DuckDB-backed) MASE, so the DuckDB query paths are exercised
# without the geometry-write path (which needs the spatial extension). Two point
# regions let us test element_type/region resolution + RI.
.lazyPointsMASE <- function() {
    skip_if_not_installed("MultiAssaySpatialExperiment")
    mat <- matrix(1:9, 3, 3,
        dimnames = list(paste0("G", 1:3), c("A", "B", "C")))
    mase <- MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
        experiments = MultiAssayExperiment::ExperimentList(assay1 = mat),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(
            assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = MultiAssaySpatialExperiment::PointsLayerList(
            centroids = S4Vectors::DataFrame(
                x = c(1.5, 2.5, 3.5), y = c(1.5, 2.5, 3.5),
                instance_id = c("A", "B", "C")),
            nuclei = S4Vectors::DataFrame(
                x = c(1.4, 2.4, 3.4), y = c(1.4, 2.4, 3.4),
                instance_id = c("A", "B", "C"))),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B", "C"),
            element_type = "points", region = "centroids",
            instance_id = c("A", "B", "C")),
        metadata = list(transforms = list(
            "points/centroids" = list(
                global = list(type = "identity"),
                scaled = list(type = "scale", scale = c(2, 2))))))
    path <- file.path(tempdir(), paste0("mase_q_", sample.int(1e6, 1L)))
    BiocDuckDB::writeParquet(mase, path)
    BiocDuckDB::readParquet(path)
}

test_that("per-element transforms round-trip and expose coordinate systems", {
    mase <- .lazyPointsMASE()
    expect_true(BiocDuckDB:::.maseHasLazyLayers(mase))
    expect_false(is.null(S4Vectors::metadata(mase)$transforms))
    expect_setequal(spatialCoordinateSystems(mase), c("global", "scaled"))
})

test_that("spatialViews registers layer + spatialMap views", {
    mase <- .lazyPointsMASE()
    v <- spatialViews(mase)
    expect_s3_class(v, "MASESpatialViews")
    expect_true("centroids" %in% names(v$points))
    expect_false(is.na(v$spatial_map))
    n <- DBI::dbGetQuery(v$conn,
        sprintf("SELECT COUNT(*) n FROM %s", v$spatial_map))[["n"]]
    expect_equal(n, 3L)
})

test_that("linkSpatialMap joins observations to their spatial coordinates", {
    mase <- .lazyPointsMASE()
    linked <- as.data.frame(linkSpatialMap(mase, assay = "assay1"))
    linked <- linked[order(linked$colname), ]
    expect_equal(nrow(linked), 3L)
    expect_equal(linked$colname, c("A", "B", "C"))
    expect_equal(linked$x, c(1.5, 2.5, 3.5))
    expect_equal(linked$y, c(1.5, 2.5, 3.5))
    # no internal __index__ leaks through
    expect_false("__index__" %in% colnames(linked))
})

test_that("linkSpatialMap requires a single (element_type, region)", {
    mase <- .lazyPointsMASE()
    # add rows referencing a second region so resolution is ambiguous
    spmap <- rbind(as.data.frame(MultiAssaySpatialExperiment::spatialMap(mase)),
        data.frame(assay = "assay1", colname = c("A", "B", "C"),
                   element_type = "points", region = "nuclei",
                   instance_id = c("A", "B", "C")))
    mase <- BiocDuckDB:::replaceSlots(mase,
        spatialMap = S4Vectors::DataFrame(spmap), check = FALSE)
    expect_error(linkSpatialMap(mase, assay = "assay1"), "single")
    # disambiguated by region -> works
    got <- as.data.frame(linkSpatialMap(mase, region = "nuclei"))
    expect_equal(sort(got$x), c(1.4, 2.4, 3.4))
})

test_that("validateSpatialMap passes a clean map and catches violations", {
    mase <- .lazyPointsMASE()
    expect_equal(nrow(validateSpatialMap(mase)), 0L)

    spmap <- rbind(as.data.frame(MultiAssaySpatialExperiment::spatialMap(mase)),
        data.frame(assay = "assay1", colname = "ZZ",  # orphan colname + instance
                   element_type = "points", region = "centroids",
                   instance_id = "ZZ"),
        data.frame(assay = "assay1", colname = "A",   # unknown layer
                   element_type = "shapes", region = "ghost",
                   instance_id = "A"))
    bad <- BiocDuckDB:::replaceSlots(mase,
        spatialMap = S4Vectors::DataFrame(spmap), check = FALSE)
    report <- validateSpatialMap(bad)
    expect_true(all(c("orphan_instance", "orphan_colname", "unknown_layer") %in%
                    report$type))
    expect_error(validateSpatialMap(bad, strict = TRUE), "referential-integrity")
})

test_that("CT-graph alignment transforms an element's coordinates (x/y path)", {
    # The #49 -> #50 wiring: build a ctGraph from metadata(mase)$transforms and
    # align a layer into a target coordinate system. Exercised on the x/y
    # arithmetic path, so no spatial extension is required.
    skip_if_not_installed("DuckDBSpatial")
    mase <- .lazyPointsMASE()
    g <- BiocDuckDB:::.maseCTgraph(mase)
    centroids <- MultiAssaySpatialExperiment::spatialPoints(mase)[["centroids"]]
    aligned <- BiocDuckDB:::.alignLayer(centroids, g, "points/centroids", "scaled")
    df <- as.data.frame(aligned)
    df <- df[order(df$x), ]
    expect_equal(df$x, c(3, 5, 7))   # 1.5, 2.5, 3.5 scaled by 2
    expect_equal(df$y, c(3, 5, 7))
    # aligning to the element's own identity target is a no-op
    same <- BiocDuckDB:::.alignLayer(centroids, g, "points/centroids", "global")
    expect_equal(sort(as.data.frame(same)$x), c(1.5, 2.5, 3.5))
})

test_that("spatialElementJoin joins across elements with CT alignment", {
    skip_if_not_installed("MultiAssaySpatialExperiment")
    skip_if_not_installed("sf")
    skip_if_not_installed("DuckDBSpatial")
    skip_if_not(tryCatch({
        DBI::dbGetQuery(DuckDBDataFrame::acquireDuckDBConn(),
            "SELECT ST_Area(ST_GeomFromText('POLYGON((0 0,0 1,1 1,1 0,0 0))'))")
        TRUE
    }, error = function(e) FALSE), "DuckDB spatial extension not available")
    mase <- makeLazySpatialMASE()  # lazy points + geometry-bearing shapes layer

    # points/centroids -> shapes/cells, point-in-polygon (same "global" frame)
    idx <- spatialElementJoin(mase, "points/centroids", "shapes/cells",
                              coordinate_system = "global")
    expect_equal(length(idx), 3L)
    expect_true("global" %in% spatialCoordinateSystems(mase))
})
