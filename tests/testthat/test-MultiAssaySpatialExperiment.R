# Tests the basic functions of a DuckDBMatrix.
# library(testthat); library(BiocDuckDB); source("setup.R"); source("test-DuckDBMatrix-spatial.R")

.make_ddb_matrix <- function(mat, label) {
    pq_dir <- file.path(tempdir(), paste0("c3_", label))
    DuckDBArray::writeCoordArray(mat, pq_dir)
    nr <- nrow(mat); nc <- ncol(mat)
    DuckDBArray::DuckDBMatrix(pq_dir, datacol = "value",
        row = list(index1 = setNames(seq_len(nr), rownames(mat))),
        col = list(index2 = setNames(seq_len(nc), colnames(mat))))
}

test_that("aggregateByRegion sum works with DuckDBMatrix assay", {
    mat <- matrix(c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9), 3, 3,
        dimnames = list(paste0("G", 1:3), c("A", "B", "C")))
    ddb_mat <- .make_ddb_matrix(mat, "sum")

    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5, 3.5),
        y = c(1.5, 2.5, 3.5),
        instance_id = c("A", "B", "C"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = c("cell1", "cell2", "cell3"),
        geometry = sf::st_sfc(
            sf::st_polygon(list(matrix(c(1,1,2,1,2,2,1,2,1,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,1,3,1,3,2,2,2,2,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,2,3,2,3,3,2,3,2,2), ncol=2, byrow=TRUE)))))

    make_mase <- function(assay_data) {
        MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
            experiments = MultiAssayExperiment::ExperimentList(assay1 = assay_data),
            colData = S4Vectors::DataFrame(row.names = "s1"),
            sampleMap = S4Vectors::DataFrame(
                assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
            points = MultiAssaySpatialExperiment::PointsLayerList(centroids = pts),
            shapes = MultiAssaySpatialExperiment::ShapesLayerList(cells = shp_df),
            spatialMap = S4Vectors::DataFrame(
                assay = "assay1", colname = c("A", "B", "C"),
                element_type = "points", region = "centroids", instance_id = c("A", "B", "C")))
    }

    se_ddb <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = ddb_mat))

    mase_mem <- MultiAssaySpatialExperiment::annotateWithRegions(make_mase(mat),
        points = "centroids", shapes = "cells")
    mase_ddb <- MultiAssaySpatialExperiment::annotateWithRegions(make_mase(se_ddb),
        points = "centroids", shapes = "cells")

    agg_mem <- MultiAssaySpatialExperiment::aggregateByRegion(mase_mem, by = "cells", FUN = "sum")
    agg_ddb <- MultiAssaySpatialExperiment::aggregateByRegion(mase_ddb, by = "cells", FUN = "sum")

    expect_equal(names(agg_ddb), names(agg_mem))
    expect_equal(agg_ddb[["assay1"]], agg_mem[["assay1"]])
})

test_that("aggregateByRegion mean works with DuckDBMatrix assay", {
    mat <- matrix(c(2.0, 4.0, 6.0, 8.0), 2, 2,
        dimnames = list(c("G1", "G2"), c("A", "B")))
    ddb_mat <- .make_ddb_matrix(mat, "mean")

    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5), y = c(1.5, 2.5), instance_id = c("A", "B"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = "cell1",
        geometry = sf::st_sfc(sf::st_polygon(list(
            matrix(c(1,1,3,1,3,3,1,3,1,1), ncol = 2, byrow = TRUE)))))

    make_mase <- function(assay_data) {
        MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
            experiments = MultiAssayExperiment::ExperimentList(assay1 = assay_data),
            colData = S4Vectors::DataFrame(row.names = "s1"),
            sampleMap = S4Vectors::DataFrame(
                assay = "assay1", primary = "s1", colname = c("A", "B")),
            points = MultiAssaySpatialExperiment::PointsLayerList(centroids = pts),
            shapes = MultiAssaySpatialExperiment::ShapesLayerList(cells = shp_df),
            spatialMap = S4Vectors::DataFrame(
                assay = "assay1", colname = c("A", "B"),
                element_type = "points", region = "centroids", instance_id = c("A", "B")))
    }

    se_ddb <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = ddb_mat))

    mase_mem <- MultiAssaySpatialExperiment::annotateWithRegions(make_mase(mat),
        points = "centroids", shapes = "cells")
    mase_ddb <- MultiAssaySpatialExperiment::annotateWithRegions(make_mase(se_ddb),
        points = "centroids", shapes = "cells")

    agg_mem <- MultiAssaySpatialExperiment::aggregateByRegion(mase_mem, by = "cells", FUN = "mean")
    agg_ddb <- MultiAssaySpatialExperiment::aggregateByRegion(mase_ddb, by = "cells", FUN = "mean")

    expect_equal(names(agg_ddb), names(agg_mem))
    expect_equal(agg_ddb[["assay1"]], agg_mem[["assay1"]])
})

test_that("aggregateByRegion count is data-agnostic (spatialMap only)", {
    mat <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), 2, 3,
        dimnames = list(c("G1", "G2"), c("A", "B", "C")))
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5), y = c(1.5, 2.5), instance_id = c("A", "B"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = "cell1",
        geometry = sf::st_sfc(sf::st_polygon(list(
            matrix(c(1,1,3,1,3,3,1,3,1,1), ncol = 2, byrow = TRUE)))))

    mase <- MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
        experiments = MultiAssayExperiment::ExperimentList(assay1 = mat),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(
            assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = MultiAssaySpatialExperiment::PointsLayerList(centroids = pts),
        shapes = MultiAssaySpatialExperiment::ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B"),
            element_type = "points", region = "centroids", instance_id = c("A", "B")))
    mase <- MultiAssaySpatialExperiment::annotateWithRegions(mase, points = "centroids", shapes = "cells")

    cnt <- MultiAssaySpatialExperiment::aggregateByRegion(mase, by = "cells", FUN = "count")
    expect_s4_class(cnt, "DataFrame")
    expect_equal(cnt$count, 2L)
})
