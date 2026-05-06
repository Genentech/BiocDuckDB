#' Spatial methods for DuckDBDataFrame
#'
#' @description
#' Methods for \code{spatialOverlaps} and \code{spatialMatch} generics defined
#' in \pkg{MultiAssaySpatialExperiment}, specialised for
#' \linkS4class{DuckDBDataFrame} objects.  These push spatial predicates
#' down to DuckDB SQL so that filtering and matching happen in the database.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{spatialOverlaps(x, y, coords = NULL, geom = "geometry", ...)}:}{
#'     When \code{coords} is provided, constructs point geometries via
#'     \code{ST_Point} and tests intersection with \code{y} in SQL.
#'     When \code{coords} is \code{NULL}, uses the geometry column named by
#'     \code{geom}.  Returns a logical vector.
#'   }
#'   \item{\code{spatialMatch(x, table, coords, geom = "geometry", join = NULL, ...)}:}{
#'     For each row in \code{x}, finds the first row in \code{table} whose
#'     geometry spatially matches.  Materialises coordinate and geometry
#'     columns and delegates to \pkg{sf} for the spatial join.  Returns an
#'     integer vector of positions in \code{table}.
#'   }
#' }
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[MultiAssaySpatialExperiment]{spatialOverlaps}} for the
#'     generic and \code{DataFrame} method.
#'   \item \code{\link[MultiAssaySpatialExperiment]{spatialMatch}} for the
#'     generic and \code{DataFrame} method.
#' }
#'
#' @aliases spatialOverlaps,DuckDBDataFrame-method
#' @aliases spatialMatch,DuckDBDataFrame,DataFrame-method
#'
#' @keywords methods
#'
#' @name DuckDBDataFrame-spatial
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Internal helper: geometry to SQL expression
###

#' @importFrom dplyr sql
.geom_to_sql <- function(y) {
    if (is.call(y))
        return(y)
    if (is.character(y))
        return(sql(sprintf("ST_GeomFromText('%s')", y)))
    if (inherits(y, "sfg"))
        return(sql(sprintf("ST_GeomFromText('%s')", sf::st_as_text(y))))
    if (inherits(y, "sfc")) {
        if (length(y) != 1L)
            stop("'y' must be a single geometry (sfc of length 1)")
        return(sql(sprintf("ST_GeomFromText('%s')", sf::st_as_text(y[[1L]]))))
    }
    stop("unsupported geometry type: ", class(y)[1L])
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### spatialOverlaps & spatialMatch
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr mutate pull sql
#' @importFrom MultiAssaySpatialExperiment spatialOverlaps
#' @importFrom stats setNames
setMethod("spatialOverlaps", "DuckDBDataFrame",
function(x, y, coords = NULL, geom = "geometry", ...) {
    y_sql <- .geom_to_sql(y)
    if (!is.null(coords)) {
        conn <- tblconn(x, select = FALSE)
        bool_sql <- sprintf( "ST_Intersects(ST_Point(\"%s\", \"%s\"), %s)",
                            coords[1L], coords[2L], as.character(y_sql))
        hit_expr <- setNames(list(sql(bool_sql)), ".hit")
        result <- pull(mutate(conn, !!!hit_expr), ".hit")
        as.logical(result)
    } else {
        geom_col <- x[[geom]]
        # Use DuckDBSpatial if available, otherwise fall back to sf
        if (requireNamespace("DuckDBSpatial", quietly = TRUE) && 
            inherits(geom_col, "DuckDBColumn")) {
            as.vector(DuckDBSpatial:::st_intersects.DuckDBColumn(geom_col, y))
        } else {
            as.vector(sf::st_intersects(geom_col, y))
        }
    }
})

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom MultiAssaySpatialExperiment spatialMatch
setMethod("spatialMatch", c("DuckDBDataFrame", "DataFrame"),
function(x, table, coords, geom = "geometry", join = NULL, ...) {
    if (is.null(join))
        join <- sf::st_intersects
    x_coords <- as.data.frame(x[, coords])
    pts_sfc <- sf::st_as_sfc(paste0("POINT(", x_coords[[coords[1L]]], " ",
                             x_coords[[coords[2L]]], ")"))
    if (is(table, "DuckDBDataFrame")) {
        tbl_geom <- as.vector(table[[geom]])
    } else {
        tbl_geom <- table[[geom]]
    }
    res <- join(pts_sfc, tbl_geom)
    n <- nrow(x)
    out <- integer(n)
    for (i in seq_len(n)) {
        idx <- res[[i]]
        out[i] <- if (length(idx) > 0L) idx[1L] else NA_integer_
    }
    out
})
