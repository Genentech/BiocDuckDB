### =========================================================================
### MultiAssaySpatialExperiment-spatial.R — lazy MASE spatial methods
### -------------------------------------------------------------------------

#' @include MultiAssaySpatialExperiments-internals.R
NULL

#' MultiAssaySpatialExperiment lazy spatial methods
#'
#' @description
#' Lazy \linkS4class{MultiAssaySpatialExperiment} methods that keep spatial
#' layers as \linkS4class{DuckDBDataFrame} objects when possible. In-memory
#' MASE defaults are used when no lazy layers are present.
#'
#' @section Parquet readers:
#' \describe{
#'   \item{\code{readParquetForMASE(file_path, col_select = NULL, ...)}:}{
#'     Files larger than 50 MB are read as lazy \code{DuckDBDataFrame} objects;
#'     smaller files use the in-memory \pkg{MultiAssaySpatialExperiment} method.
#'   }
#'   \item{\code{readGeoParquetForMASE(file_path, ...)}:}{
#'     Uses \code{\link[DuckDBSpatial]{readGeoParquet}} when \pkg{DuckDBSpatial}
#'     is installed; otherwise falls back to \pkg{sfarrow}.
#'   }
#' }
#'
#' @section Spatial operations:
#' \describe{
#'   \item{\code{annotateWithRegions(x, points, shapes, ...)}:}{
#'     Annotate points with shape regions using lazy \code{spatialMatch}.
#'   }
#'   \item{\code{subsetByBoundingBox(x, xmin, xmax, ymin, ymax, ...)}:}{
#'     Subset by bounding box with lazy layer filtering when applicable.
#'   }
#'   \item{\code{subsetByPolygon(x, polygon, ...)}:}{
#'     Subset by polygon with lazy layer filtering when applicable.
#'   }
#' }
#'
#' @seealso
#' \code{\link[MultiAssaySpatialExperiment]{readParquetForMASE}},
#' \code{\link[MultiAssaySpatialExperiment]{readGeoParquetForMASE}},
#' \code{\link[MultiAssaySpatialExperiment]{annotateWithRegions}},
#' \code{\link[MultiAssaySpatialExperiment]{subsetByBoundingBox}}, and
#' \code{\link[MultiAssaySpatialExperiment]{subsetByPolygon}}.
#'
#' @aliases readParquetForMASE,character-method
#' @aliases readGeoParquetForMASE,character-method
#' @aliases annotateWithRegions,MultiAssaySpatialExperiment-method
#' @aliases subsetByBoundingBox,MultiAssaySpatialExperiment-method
#' @aliases subsetByPolygon,MultiAssaySpatialExperiment-method
#'
#' @name MultiAssaySpatialExperiment-spatial
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### readParquetForMASE / readGeoParquetForMASE
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom MultiAssaySpatialExperiment readParquetForMASE readGeoParquetForMASE
#' @importFrom S4Vectors wmsg

.MASE_readParquetForMASE <- .captureMASEMethod("readParquetForMASE", "character")

#' @export
#' @importFrom DuckDBDataFrame DuckDBDataFrame
setMethod("readParquetForMASE", "character",
function(file_path, col_select = NULL, ...) {
    size_mb <- if (file.exists(file_path)) file.info(file_path)[["size"]] / (1024 * 1024) else 0
    if (size_mb > 50) {
        DuckDBDataFrame(file_path)
    } else {
        .MASE_readParquetForMASE(file_path, col_select = col_select, ...)
    }
})

#' @export
setMethod("readGeoParquetForMASE", "character",
function(file_path, ...) {
    if (requireNamespace("DuckDBSpatial", quietly = TRUE)) {
        DuckDBSpatial::readGeoParquet(file_path, ...)
    } else {
        if (!requireNamespace("sfarrow", quietly = TRUE))
            stop(wmsg("Package 'DuckDBSpatial' or 'sfarrow' required for GeoParquet reading"))
        sfarrow::st_read_parquet(file_path, ...)
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### annotateWithRegions
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom MultiAssaySpatialExperiment MultiAssaySpatialExperiment
#' @importFrom MultiAssaySpatialExperiment annotateWithRegions spatialMatch
#' @importFrom MultiAssaySpatialExperiment spatialPoints spatialShapes spatialMap

.MASE_annotateWithRegions <- .captureMASEMethod("annotateWithRegions")

#' @export
setMethod("annotateWithRegions", "MultiAssaySpatialExperiment",
function(x, points, shapes, spatialCoordsNames = c("x", "y"), regionCol = NULL,
    join = NULL) {
    if (!.maseHasLazyLayers(x))
        return(.MASE_annotateWithRegions(x, points, shapes,
            spatialCoordsNames = spatialCoordsNames, regionCol = regionCol, join = join))
    if (is.null(join))
        join <- sf::st_intersects
    pts_el <- spatialPoints(x)[[points]]
    shps_el <- spatialShapes(x)[[shapes]]
    spmap <- spatialMap(x)
    if (is.null(pts_el) || is.null(shps_el) || is.null(spmap) || nrow(spmap) == 0L)
        return(x)
    if (is.null(regionCol))
        regionCol <- shapes
    x_col <- spatialCoordsNames[1L]
    y_col <- spatialCoordsNames[2L]
    idx <- spatialMatch(pts_el, shps_el, coords = c(x_col, y_col),
                        geom = "geometry", join = join)
    shp_ids <- shps_el[["instance_id"]]
    mapped <- ifelse(is.na(idx), NA, as.vector(shp_ids)[idx])
    keep <- spmap[["region"]] == points
    if (!any(keep))
        return(x)
    new_col <- spmap[[regionCol]]
    if (is.null(new_col))
        new_col <- rep(NA, nrow(spmap))
    pt_ids <- as.character(spmap[["instance_id"]][keep])
    if ("instance_id" %in% colnames(pts_el)) {
        pt_map <- setNames(mapped, as.character(as.vector(pts_el[["instance_id"]])))
        new_col[keep] <- pt_map[pt_ids]
    } else {
        new_col[keep] <- mapped[match(pt_ids, seq_along(mapped))]
    }
    spmap[[regionCol]] <- new_col
    replaceSlots(x, spatialMap = spmap, check = FALSE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subsetByBoundingBox / subsetByPolygon
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importClassesFrom MultiAssaySpatialExperiment MultiAssaySpatialExperiment
#' @importFrom MultiAssaySpatialExperiment subsetByBoundingBox subsetByPolygon
#' @importFrom MultiAssaySpatialExperiment spatialPoints spatialShapes
#' @importFrom MultiAssaySpatialExperiment PointsLayerList ShapesLayerList
#' @importFrom MultiAssaySpatialExperiment spatialOverlaps

.MASE_subsetByBoundingBox <- .captureMASEMethod("subsetByBoundingBox")
.MASE_subsetByPolygon <- .captureMASEMethod("subsetByPolygon")

.subsetLayerByBbox <- function(el, xmin, xmax, ymin, ymax, x_col, y_col) {
    if (.isLazySpatialLayer(el)) {
        .requireDuckDBSpatial("lazy spatial subset")
        DuckDBSpatial::layerSubsetByBbox(el, xmin, xmax, ymin, ymax,
                                         x_col = x_col, y_col = y_col)
    } else {
        xv <- el[[x_col]]; yv <- el[[y_col]]
        which((xv >= xmin & xv <= xmax) & (yv >= ymin & yv <= ymax))
    }
}

.subsetLayerByBboxShapes <- function(el, xmin, xmax, ymin, ymax) {
    if (!"geometry" %in% colnames(el))
        return(integer(0L))
    if (.isLazySpatialLayer(el)) {
        .requireDuckDBSpatial("lazy spatial subset")
        DuckDBSpatial::layerSubsetByBbox(el, xmin, xmax, ymin, ymax,
                                         x_col = NULL, geom = "geometry")
    } else {
        env <- sf::st_as_sfc(sf::st_bbox(c(xmin = xmin, ymin = ymin,
                                            xmax = xmax, ymax = ymax)))
        which(spatialOverlaps(el, env, geom = "geometry"))
    }
}

.subsetLayerByGeometry <- function(el, y, coords = NULL, geom = "geometry") {
    if (.isLazySpatialLayer(el)) {
        .requireDuckDBSpatial("lazy spatial subset")
        DuckDBSpatial::layerSubsetByGeometry(el, y, coords = coords, geom = geom)
    } else {
        which(spatialOverlaps(el, y, coords = coords, geom = geom))
    }
}

.getInstanceIdsFromLayer <- function(el, idx) {
    if (length(idx) == 0L)
        return(character(0L))
    if ("instance_id" %in% colnames(el))
        as.character(el[["instance_id"]][idx])
    else
        as.character(idx)
}

#' @export
setMethod("subsetByBoundingBox", "MultiAssaySpatialExperiment",
function(x, xmin, xmax, ymin, ymax, x_col = "x", y_col = "y", ...) {
    if (!.maseHasLazyLayers(x))
        return(.MASE_subsetByBoundingBox(x, xmin, xmax, ymin, ymax,
            x_col = x_col, y_col = y_col, ...))
    pts <- spatialPoints(x)
    shps <- spatialShapes(x)
    retained <- list()
    for (nm in names(pts)) {
        el <- pts[[nm]]
        idx <- .subsetLayerByBbox(el, xmin, xmax, ymin, ymax, x_col, y_col)
        retained[[nm]] <- .getInstanceIdsFromLayer(el, idx)
    }
    for (nm in names(shps)) {
        el <- shps[[nm]]
        idx <- .subsetLayerByBboxShapes(el, xmin, xmax, ymin, ymax)
        if (length(idx))
            retained[[nm]] <- .getInstanceIdsFromLayer(el, idx)
    }
    x <- .propagateSpatialSubset(x, retained)
    filt_pts <- lapply(names(pts), function(nm) {
        el <- pts[[nm]]
        el[.subsetLayerByBbox(el, xmin, xmax, ymin, ymax, x_col, y_col), , drop = FALSE]
    })
    names(filt_pts) <- names(pts)
    filt_shps <- lapply(names(shps), function(nm) {
        el <- shps[[nm]]
        el[.subsetLayerByBboxShapes(el, xmin, xmax, ymin, ymax), , drop = FALSE]
    })
    names(filt_shps) <- names(shps)
    replaceSlots(x,
        points = PointsLayerList(filt_pts),
        shapes = ShapesLayerList(filt_shps),
        check = FALSE)
})

#' @export
setMethod("subsetByPolygon", "MultiAssaySpatialExperiment",
function(x, polygon, x_col = "x", y_col = "y", ...) {
    if (!.maseHasLazyLayers(x))
        return(.MASE_subsetByPolygon(x, polygon, x_col = x_col, y_col = y_col, ...))
    if (inherits(polygon, "sfg"))
        polygon <- sf::st_sfc(polygon)
    else if (inherits(polygon, "sf"))
        polygon <- sf::st_geometry(polygon)
    pts <- spatialPoints(x)
    shps <- spatialShapes(x)
    retained <- list()
    for (nm in names(pts)) {
        el <- pts[[nm]]
        idx <- .subsetLayerByGeometry(el, polygon, coords = c(x_col, y_col))
        retained[[nm]] <- .getInstanceIdsFromLayer(el, idx)
    }
    for (nm in names(shps)) {
        el <- shps[[nm]]
        idx <- .subsetLayerByGeometry(el, polygon, geom = "geometry")
        if (length(idx))
            retained[[nm]] <- .getInstanceIdsFromLayer(el, idx)
    }
    x <- .propagateSpatialSubset(x, retained)
    filt_pts <- lapply(names(pts), function(nm) {
        el <- pts[[nm]]
        el[.subsetLayerByGeometry(el, polygon, coords = c(x_col, y_col)), , drop = FALSE]
    })
    names(filt_pts) <- names(pts)
    filt_shps <- lapply(names(shps), function(nm) {
        el <- shps[[nm]]
        el[.subsetLayerByGeometry(el, polygon, geom = "geometry"), , drop = FALSE]
    })
    names(filt_shps) <- names(shps)
    replaceSlots(x,
        points = PointsLayerList(filt_pts),
        shapes = ShapesLayerList(filt_shps),
        check = FALSE)
})
