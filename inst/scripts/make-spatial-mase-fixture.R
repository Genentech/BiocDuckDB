## Generate the minimal spatial MASE fixture shipped in inst/extdata/spatial_mase/.
## It is used by the runnable examples of the spatial-map query functions
## (linkSpatialMap, spatialCoordinateSystems, spatialViews, validateSpatialMap),
## so those examples are self-contained and runnable at R CMD check time without
## a large real dataset.
##
## The fixture is a points-only, DuckDB-backed MultiAssaySpatialExperiment with a
## spatial map linking observations to point coordinates, plus two per-element
## coordinate systems (an identity "global" and a scaled one), serialized with
## writeParquet().
##
## Regenerate from the package root with:
##   Rscript inst/scripts/make-spatial-mase-fixture.R

suppressMessages(library(BiocDuckDB))

mat <- matrix(1:9, 3, 3, dimnames = list(paste0("G", 1:3), c("A", "B", "C")))
mase <- MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
    experiments = MultiAssayExperiment::ExperimentList(assay1 = mat),
    colData = S4Vectors::DataFrame(row.names = "s1"),
    sampleMap = S4Vectors::DataFrame(assay = "assay1", primary = "s1",
        colname = c("A", "B", "C")),
    points = MultiAssaySpatialExperiment::PointsLayerList(
        centroids = S4Vectors::DataFrame(x = c(1.5, 2.5, 3.5),
            y = c(1.5, 2.5, 3.5), instance_id = c("A", "B", "C")),
        nuclei = S4Vectors::DataFrame(x = c(1.4, 2.4, 3.4),
            y = c(1.4, 2.4, 3.4), instance_id = c("A", "B", "C"))),
    spatialMap = S4Vectors::DataFrame(assay = "assay1",
        colname = c("A", "B", "C"), element_type = "points",
        region = "centroids", instance_id = c("A", "B", "C")),
    metadata = list(transforms = list("points/centroids" = list(
        global = list(type = "identity"),
        scaled = list(type = "scale", scale = c(2, 2))))))

out <- file.path("inst", "extdata", "spatial_mase")
unlink(out, recursive = TRUE)
writeParquet(mase, out)
