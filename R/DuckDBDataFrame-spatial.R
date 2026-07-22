### =========================================================================
### DuckDBDataFrame-spatial.R — MASE generic registration for lazy layers
### -------------------------------------------------------------------------
### Spatial engines live in DuckDBSpatial; BiocDuckDB registers MASE methods.

#' @include MultiAssaySpatialExperiments-internals.R
NULL

#' DuckDBDataFrame spatial methods
#'
#' @description
#' \linkS4class{DuckDBDataFrame} methods for \pkg{MultiAssaySpatialExperiment}
#' spatial generics. Lazy layers delegate to \pkg{DuckDBSpatial} for
#' SQL-backed spatial operations.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{spatialOverlaps(x, y, coords = NULL, geom = "geometry")}:}{
#'     Test spatial overlap between lazy spatial layers via
#'     \code{\link[DuckDBSpatial]{layerSpatialOverlaps}}.
#'   }
#'   \item{\code{spatialMatch(x, table, coords, geom = "geometry", join = NULL)}:}{
#'     Match points or geometries to a table via
#'     \code{\link[DuckDBSpatial]{layerSpatialMatch}}.
#'   }
#'   \item{\code{spatialJoin(x, y, join = st_intersects, sparse = TRUE)}:}{
#'     Spatial join between lazy layers. When \code{sparse = TRUE}, returns a
#'     sparse pair table via \code{\link[DuckDBSpatial]{st_intersects_table}};
#'     otherwise uses \code{\link[sf:st_join]{st_join}}.
#'   }
#' }
#'
#' @seealso
#' \code{\link[MultiAssaySpatialExperiment]{spatialOverlaps}},
#' \code{\link[MultiAssaySpatialExperiment]{spatialMatch}}, and
#' \code{\link[MultiAssaySpatialExperiment]{spatialJoin}} for generic
#' documentation.
#'
#' @examples
#' # Spatial predicates need DuckDB's optional spatial extension. Probe for it so
#' # the example runs only where the extension can actually load; it is skipped
#' # cleanly when offline or behind a firewall that blocks the extension
#' # repository, or when DuckDBSpatial is not installed.
#' spatial_ok <- requireNamespace("DuckDBSpatial", quietly = TRUE) &&
#'     tryCatch({
#'         DBI::dbGetQuery(DuckDBDataFrame::acquireDuckDBConn(),
#'             "SELECT ST_Area(ST_GeomFromText('POLYGON((0 0,0 1,1 1,1 0,0 0))'))")
#'         TRUE
#'     }, error = function(e) FALSE)
#' if (spatial_ok) {
#'     df <- data.frame(id = 1:3, x = 1:3, y = 1:3)
#'     path <- tempfile(fileext = ".parquet")
#'     arrow::write_parquet(df, path)
#'     ddb <- DuckDBDataFrame(path, datacols = c("x", "y"), keycol = "id")
#'     polygon <- "POLYGON((0 0, 6 0, 6 6, 0 6, 0 0))"
#'     spatialOverlaps(ddb, polygon, coords = c("x", "y"))
#'     unlink(path)
#' }
#'
#' @return
#' \code{spatialOverlaps()} returns the pairwise spatial-overlap result between
#' the lazy layers, \code{spatialMatch()} the matched table, and
#' \code{spatialJoin()} a sparse pair table (mapping row indices between layers)
#' when \code{sparse = TRUE}, or the joined table otherwise. All are computed
#' lazily over the DuckDB-backed layers.
#'
#' @aliases spatialOverlaps,DuckDBDataFrame-method
#' @aliases spatialMatch,DuckDBDataFrame,DataFrame-method
#' @aliases spatialMatch,DuckDBDataFrame,DuckDBDataFrame-method
#' @aliases spatialJoin,DuckDBDataFrame,DuckDBDataFrame-method
#'
#' @name DuckDBDataFrame-spatial
NULL

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom MultiAssaySpatialExperiment spatialOverlaps
setMethod("spatialOverlaps", "DuckDBDataFrame",
function(x, y, coords = NULL, geom = "geometry", ...) {
    .requireDuckDBSpatial("spatialOverlaps on DuckDBDataFrame")
    DuckDBSpatial::layerSpatialOverlaps(x, y, coords = coords, geom = geom)
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom MultiAssaySpatialExperiment spatialMatch
setMethod("spatialMatch", c("DuckDBDataFrame", "DataFrame"),
function(x, table, coords, geom = "geometry", join = NULL, ...) {
    .requireDuckDBSpatial("spatialMatch on DuckDBDataFrame")
    DuckDBSpatial::layerSpatialMatch(x, table, coords = coords, geom = geom, join = join)
})

#' @export
setMethod("spatialMatch", c("DuckDBDataFrame", "DuckDBDataFrame"),
function(x, table, coords, geom = "geometry", join = NULL, ...) {
    .requireDuckDBSpatial("spatialMatch on DuckDBDataFrame")
    DuckDBSpatial::layerSpatialMatch(x, table, coords = coords, geom = geom, join = join)
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom MultiAssaySpatialExperiment spatialJoin
setMethod("spatialJoin", c("DuckDBDataFrame", "DuckDBDataFrame"),
function(x, y, join = sf::st_intersects, sparse = TRUE, ...) {
    .requireDuckDBSpatial("spatialJoin on DuckDBDataFrame")
    if (isTRUE(sparse))
        DuckDBSpatial::st_intersects_table(x, y)
    else
        sf::st_join(x, y, join = join, ...)
})
