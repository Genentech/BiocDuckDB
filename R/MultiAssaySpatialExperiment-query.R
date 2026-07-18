### =========================================================================
### MultiAssaySpatialExperiment-query.R — cross-element DuckDB query layer
### -------------------------------------------------------------------------
###
### A relational DuckDB query layer over a MultiAssaySpatialExperiment (MASE)
### whose spatial layers are DuckDB-backed (DuckDBDataFrame). It exposes each
### element and the spatialMap junction as on-the-fly DuckDB views, links assay
### observations to spatial layer rows through the full spatialMap key
### (assay, colname, element_type, region, instance_id), validates referential
### integrity with DuckDB anti-joins, and runs ST_* cross-element joins (pushed
### to DuckDBSpatial) after aligning elements into a common coordinate system
### via the coordinate-transform graph (DuckDBSpatial, ROADMAP #49). MASE itself
### stays DuckDB-free; this is BiocDuckDB's DuckDB engine operating on a MASE.
###
### Per-element transforms live in metadata(mase)$transforms as
###   "<element_type>/<region>" -> list(<coordinate_system> -> <transform dict>)
### (RFC-5 / scibis-shaped dicts, e.g. list(type = "scale", scale = c(2, 3))),
### and round-trip through the standard MASE metadata annotations.

#' @include MultiAssaySpatialExperiments-internals.R
NULL

# Unique temp relation name (R-side randomness; mirrors
# writeSpatialPointsParquet).
#' @importFrom stats runif
.maseTmpName <- function(prefix) {
    sprintf("%s%d", prefix, as.integer(runif(1L, 1L, .Machine$integer.max)))
}

.maseTransforms <- function(mase) {
    tf <- S4Vectors::metadata(mase)[["transforms"]]
    if (is.null(tf)) list() else tf
}

#' Coordinate systems of a MultiAssaySpatialExperiment
#'
#' Returns the names of the coordinate systems referenced by the per-element
#' transforms stored in \code{metadata(x)$transforms} (the targets each element
#' is mapped into). Empty when no transforms are recorded.
#'
#' @param x A \linkS4class{MultiAssaySpatialExperiment}.
#' @return A character vector of coordinate-system names.
#' @examples
#' # metadata(mase)$transforms <- list("points/centroids" =
#' #     list(global = list(type = "identity")))
#' # spatialCoordinateSystems(mase)  # -> "global"
#' @export
#' @importFrom S4Vectors metadata
spatialCoordinateSystems <- function(x) {
    tf <- .maseTransforms(x)
    if (!length(tf)) return(character(0L))
    unique(unlist(lapply(tf, names), use.names = FALSE))
}

# Build a DuckDBSpatial coordinate-transform graph from a MASE's transforms.
# Element intrinsic coordinate systems are the "<type>/<region>" node names;
# each stored transform is an edge element -> coordinate_system. All nodes are
# 2-D (x, y) -- the DuckDB apply-path is 2-D.
.maseCTgraph <- function(mase) {
    .requireDuckDBSpatial("coordinate-transform alignment")
    tf <- .maseTransforms(mase)
    if (!length(tf))
        stop("no coordinate transforms in metadata(mase)$transforms")
    edges <- list()
    for (elem in names(tf)) {
        for (cs in names(tf[[elem]])) {
            edges <- c(edges, list(list(
                input = elem, output = cs,
                transform = DuckDBSpatial::asCoordinateTransform(tf[[elem]][[cs]]))))
        }
    }
    nodes <- unique(c(vapply(edges, `[[`, character(1L), "input"),
                      vapply(edges, `[[`, character(1L), "output")))
    systems <- stats::setNames(rep(list(c("x", "y")), length(nodes)), nodes)
    DuckDBSpatial::ctGraph(edges, systems = systems)
}

.elementKey <- function(element_type, region) paste0(element_type, "/", region)

# Resolve a "<type>/<region>" (or list(element_type=, region=)) spec to its
# layer.
.resolveElement <- function(mase, spec) {
    if (is.list(spec)) {
        et <- spec[["element_type"]]; rg <- spec[["region"]]
    } else {
        parts <- strsplit(spec, "/", fixed = TRUE)[[1L]]
        if (length(parts) != 2L)
            stop("element spec must be '<element_type>/<region>', got '", spec,
                 "'")
        et <- parts[1L]; rg <- parts[2L]
    }
    layer <- switch(et,
        points = MultiAssaySpatialExperiment::spatialPoints(mase)[[rg]],
        shapes = MultiAssaySpatialExperiment::spatialShapes(mase)[[rg]],
        stop("element_type must be 'points' or 'shapes', got '", et, "'"))
    if (is.null(layer))
        stop("no ", et, " layer named '", rg, "'")
    list(layer = layer, element_type = et, region = rg,
         key = .elementKey(et, rg))
}

# Register one layer as a DuckDB temp view (lazy layer -> view over its rendered
# SQL; materialized layer -> temp table) and return the view name.
#' @importFrom DBI dbExecute dbWriteTable dbQuoteIdentifier
#' @importFrom dbplyr sql_render
.registerLayerView <- function(layer, conn, vname) {
    if (.isLazySpatialLayer(layer)) {
        sql <- as.character(sql_render(DuckDBDataFrame::tblconn(layer, select = FALSE)))
        dbExecute(conn, sprintf("CREATE OR REPLACE TEMP VIEW %s AS %s",
                                dbQuoteIdentifier(conn, vname), sql))
    } else {
        dbWriteTable(conn, vname, as.data.frame(layer),
                     temporary = TRUE, overwrite = TRUE)
    }
    vname
}

#' On-the-fly DuckDB views over a MASE's spatial elements
#'
#' Registers each spatial layer (points, shapes) and the \code{spatialMap}
#' junction of a \linkS4class{MultiAssaySpatialExperiment} as DuckDB temp views
#' on a shared connection, so cross-element queries can be expressed as SQL.
#' Lazy (\linkS4class{DuckDBDataFrame}) layers are registered as views over
#' their rendered SQL (no materialization); in-memory layers and the
#' materialized \code{spatialMap} are registered as temp tables. This is the
#' substrate used by \code{\link{linkSpatialMap}} /
#' \code{\link{validateSpatialMap}} and a handle for power-user raw SQL.
#'
#' @param mase A \linkS4class{MultiAssaySpatialExperiment}.
#' @param conn A DuckDB connection (default the shared BiocDuckDB connection).
#' @param prefix View-name prefix.
#' @return A \code{MASESpatialViews} registry (a list of view names + \code{conn}).
#' @examples
#' # v <- spatialViews(mase)
#' # DBI::dbGetQuery(v$conn, sprintf("SELECT * FROM %s LIMIT 5", v$spatial_map))
#' @export
#' @importFrom DuckDBDataFrame acquireDuckDBConn
#' @importFrom DBI dbWriteTable
spatialViews <- function(mase, conn = acquireDuckDBConn(), prefix = "mase_") {
    reg <- list(conn = conn, points = character(0L), shapes = character(0L),
                spatial_map = NA_character_)
    for (slot in c("points", "shapes")) {
        layers <- if (slot == "points")
            MultiAssaySpatialExperiment::spatialPoints(mase)
        else MultiAssaySpatialExperiment::spatialShapes(mase)
        vnames <- character(0L)
        for (nm in names(layers)) {
            el <- layers[[nm]]
            if (is.null(el) || !nrow(el)) next
            vnames[[nm]] <- .registerLayerView(el, conn,
                sprintf("%s%s_%s", prefix, slot, nm))
        }
        reg[[slot]] <- vnames
    }
    spmap <- MultiAssaySpatialExperiment::spatialMap(mase)
    if (!is.null(spmap) && nrow(spmap)) {
        vname <- paste0(prefix, "spatial_map")
        dbWriteTable(conn, vname, as.data.frame(spmap),
                     temporary = TRUE, overwrite = TRUE)
        reg[["spatial_map"]] <- vname
    }
    structure(reg, class = "MASESpatialViews")
}

#' @export
print.MASESpatialViews <- function(x, ...) {
    cat("<MASESpatialViews: ", length(x$points), " points, ", length(x$shapes),
        " shapes, spatial_map=",
        if (is.na(x$spatial_map)) "none" else x$spatial_map, ">\n", sep = "")
    invisible(x)
}

#' Link assay observations to their spatial layer rows
#'
#' Joins the observations recorded in \code{spatialMap} to the rows of the
#' spatial layer they reference, through the full junction key
#' (\code{assay}, \code{colname}, \code{element_type}, \code{region},
#' \code{instance_id}). Resolves a single \code{(element_type, region)} layer
#' (inferred from the filtered \code{spatialMap} when unambiguous, else
#' specify); the join is pushed to DuckDB and returned as a lazy
#' \linkS4class{DuckDBDataFrame} with the assay identity columns plus the
#' layer's coordinate columns.
#'
#' @param mase A \linkS4class{MultiAssaySpatialExperiment}.
#' @param assay,element_type,region Optional filters selecting the spatialMap
#'   rows (and thus the single layer) to link.
#' @param conn A DuckDB connection (default the shared BiocDuckDB connection).
#' @return A lazy \linkS4class{DuckDBDataFrame} of linked observations.
#' @examples
#' # linkSpatialMap(mase, assay = "assay1")
#' @export
#' @importFrom DuckDBDataFrame acquireDuckDBConn DuckDBDataFrame
#' @importFrom DBI dbExecute dbWriteTable dbQuoteIdentifier
linkSpatialMap <-
function(mase, assay = NULL, element_type = NULL, region = NULL,
         conn = acquireDuckDBConn())
{
    spmap <- as.data.frame(MultiAssaySpatialExperiment::spatialMap(mase))
    if (is.null(spmap) || !nrow(spmap))
        stop("'mase' has an empty spatialMap")
    if (!is.null(assay)) spmap <- spmap[spmap[["assay"]] == assay, , drop = FALSE]
    if (!is.null(element_type))
        spmap <- spmap[spmap[["element_type"]] == element_type, , drop = FALSE]
    if (!is.null(region)) spmap <- spmap[spmap[["region"]] == region, , drop = FALSE]
    if (!nrow(spmap))
        stop("no spatialMap rows match the given assay/element_type/region")
    combos <- unique(spmap[c("element_type", "region")])
    if (nrow(combos) > 1L)
        stop("linkSpatialMap resolves a single (element_type, region); specify ",
             "element_type/region. Present: ",
             paste(sprintf("%s/%s", combos$element_type, combos$region),
                   collapse = ", "))
    el <- .resolveElement(mase, list(element_type = combos$element_type[1L],
                                     region = combos$region[1L]))
    lv <- .registerLayerView(el$layer, conn, .maseTmpName("mase_link_layer_"))
    smv <- .maseTmpName("mase_link_spmap_")
    dbWriteTable(conn, smv, spmap, temporary = TRUE, overwrite = TRUE)
    # Select the layer's declared data columns (x/y/geometry, ...), not the
    # internal __index__ / keycol columns carried in the rendered relation.
    layer_cols <- setdiff(colnames(el$layer), "instance_id")
    sel_layer <- paste(sprintf("c.%s", dbQuoteIdentifier(conn, layer_cols)),
                       collapse = ", ")
    join_view <- .maseTmpName("mase_link_")
    sql <- sprintf(paste0(
        "CREATE OR REPLACE TEMP VIEW %s AS ",
        "SELECT m.assay, m.colname, m.element_type, m.region, m.instance_id, %s ",
        "FROM %s m JOIN %s c ",
        "ON CAST(m.instance_id AS VARCHAR) = CAST(c.instance_id AS VARCHAR)"),
        dbQuoteIdentifier(conn, join_view), sel_layer,
        dbQuoteIdentifier(conn, smv), dbQuoteIdentifier(conn, lv))
    dbExecute(conn, sql)
    DuckDBDataFrame(dplyr::tbl(conn, join_view))
}

#' Validate a MASE's spatialMap referential integrity
#'
#' Checks the \code{spatialMap} foreign keys against the spatial layers and the
#' assays, out-of-core via DuckDB anti-joins on lazy layers. Reports four
#' violation classes: \code{unknown_layer} (a \code{(element_type, region)} with
#' no matching layer), \code{orphan_instance} (a \code{spatialMap} row whose
#' \code{instance_id} is absent from its layer), \code{orphan_colname} (an
#' \code{(assay, colname)} that is not a column of that experiment), and
#' \code{duplicate_instance} (an \code{instance_id} occurring more than once
#' within a layer).
#'
#' @param mase A \linkS4class{MultiAssaySpatialExperiment}.
#' @param strict If \code{TRUE}, error when any violation is found instead of
#'   returning the report.
#' @param conn A DuckDB connection (default the shared BiocDuckDB connection).
#' @return A \code{DataFrame} of violations (\code{type}, \code{assay},
#'   \code{colname}, \code{element_type}, \code{region}, \code{instance_id},
#'   \code{detail}); empty when the spatialMap is valid.
#' @examples
#' # validateSpatialMap(mase)              # report
#' # validateSpatialMap(mase, strict = TRUE)  # error on any violation
#' @export
#' @importFrom DuckDBDataFrame acquireDuckDBConn
#' @importFrom DBI dbExecute dbWriteTable dbGetQuery dbQuoteIdentifier
#' @importFrom S4Vectors DataFrame
validateSpatialMap <-
function(mase, strict = FALSE, conn = acquireDuckDBConn())
{
    empty <- DataFrame(type = character(0L), assay = character(0L),
                       colname = character(0L), element_type = character(0L),
                       region = character(0L), instance_id = character(0L),
                       detail = character(0L))
    spmap <- MultiAssaySpatialExperiment::spatialMap(mase)
    if (is.null(spmap) || !nrow(spmap)) return(empty)
    spmap <- as.data.frame(spmap)
    viol <- list()
    add <- function(rows, type, detail = NA_character_) {
        if (!nrow(rows)) return()
        viol[[length(viol) + 1L]] <<- DataFrame(
            type = type, assay = as.character(rows[["assay"]]),
            colname = as.character(rows[["colname"]]),
            element_type = as.character(rows[["element_type"]]),
            region = as.character(rows[["region"]]),
            instance_id = as.character(rows[["instance_id"]]), detail = detail)
    }

    # (3) orphan (assay, colname): not a column of that experiment (in-memory)
    exps <- MultiAssayExperiment::experiments(mase)
    for (a in unique(spmap[["assay"]])) {
        rows <- spmap[spmap[["assay"]] == a, , drop = FALSE]
        cols <- if (a %in% names(exps)) colnames(exps[[a]]) else character(0L)
        bad <- rows[!(as.character(rows[["colname"]]) %in% cols), , drop = FALSE]
        add(bad, "orphan_colname",
            if (a %in% names(exps)) "colname not in experiment" else "unknown assay")
    }

    # per (element_type, region): unknown layer / orphan + duplicate instance_id
    combos <- unique(spmap[c("element_type", "region")])
    for (i in seq_len(nrow(combos))) {
        et <- combos$element_type[i]; rg <- combos$region[i]
        rows <- spmap[spmap[["element_type"]] == et & spmap[["region"]] == rg, ,
                      drop = FALSE]
        layer <- tryCatch(.resolveElement(mase, list(element_type = et,
                                                     region = rg))$layer,
                          error = function(e) NULL)
        if (is.null(layer)) {
            add(rows, "unknown_layer", "no layer for element_type/region")
            next
        }
        lv <- .registerLayerView(layer, conn, .maseTmpName("mase_ri_layer_"))
        smv <- .maseTmpName("mase_ri_spmap_")
        dbWriteTable(conn, smv, rows, temporary = TRUE, overwrite = TRUE)
        orphan <- dbGetQuery(conn, sprintf(paste0(
            "SELECT m.* FROM %s m LEFT JOIN %s c ",
            "ON CAST(m.instance_id AS VARCHAR) = CAST(c.instance_id AS VARCHAR) ",
            "WHERE c.instance_id IS NULL"),
            dbQuoteIdentifier(conn, smv), dbQuoteIdentifier(conn, lv)))
        add(orphan, "orphan_instance", "instance_id not in layer")
        if ("instance_id" %in% colnames(layer)) {
            dups <- dbGetQuery(conn, sprintf(paste0(
                "SELECT CAST(instance_id AS VARCHAR) AS instance_id, COUNT(*) n ",
                "FROM %s GROUP BY 1 HAVING COUNT(*) > 1"),
                dbQuoteIdentifier(conn, lv)))
            if (nrow(dups)) {
                add(data.frame(assay = NA_character_, colname = NA_character_,
                               element_type = et, region = rg,
                               instance_id = dups[["instance_id"]]),
                    "duplicate_instance", "instance_id repeated in layer")
            }
        }
    }

    report <- if (length(viol)) do.call(rbind, viol) else empty
    if (isTRUE(strict) && nrow(report))
        stop("spatialMap failed referential-integrity validation: ",
             nrow(report), " violation(s); see validateSpatialMap(strict=FALSE)")
    report
}

# Align a layer into a target coordinate system via the CT graph (#49).
.alignLayer <- function(layer, graph, from, to, x_col = "x", y_col = "y",
                        geom = "geometry") {
    tf <- DuckDBSpatial::ctPath(graph, from, to)
    if (identical(tf$type, "identity")) return(layer)
    has_geom <- geom %in% colnames(layer)
    DuckDBSpatial::transformLayer(layer, tf,
        x_col = if (has_geom) NULL else x_col,
        y_col = if (has_geom) NULL else y_col, geom = geom)
}

#' Spatial join between two elements of a MASE
#'
#' Joins two spatial elements (points or shapes, each named
#' \code{"<element_type>/<region>"}) by a spatial predicate, pushed to DuckDB
#' via \pkg{DuckDBSpatial}. When \code{coordinate_system} is given and the MASE
#' carries transforms, each element is first aligned into that coordinate system
#' through the coordinate-transform graph
#' (\code{\link{spatialCoordinateSystems}}), so elements defined in different
#' intrinsic frames are compared correctly. This generalizes
#' \code{annotateWithRegions} (point-in-polygon only) to arbitrary
#' element pairs and predicates.
#'
#' @param mase A \linkS4class{MultiAssaySpatialExperiment}.
#' @param x,y Element specs \code{"<element_type>/<region>"} (e.g.
#'   \code{"points/centroids"}), or \code{list(element_type=, region=)}.
#' @param join A spatial predicate function (default
#'   \code{\link[sf:geos_binary_pred]{st_intersects}}).
#' @param coordinate_system Optional target coordinate-system name to align both
#'   elements into before joining.
#' @param x_col,y_col Coordinate column names for point elements.
#' @return The spatial-join result from \pkg{DuckDBSpatial}
#'   (\code{\link{spatialMatch}} indices for a point/geometry query).
#' @examples
#' # spatialElementJoin(mase, "points/centroids", "shapes/cells")
#' # spatialElementJoin(mase, "points/centroids", "shapes/cells",
#' #     coordinate_system = "global")
#' @export
#' @importFrom MultiAssaySpatialExperiment spatialMatch
spatialElementJoin <- function(mase, x, y, join = sf::st_intersects,
                               coordinate_system = NULL, x_col = "x",
                               y_col = "y") {
    .requireDuckDBSpatial("spatialElementJoin")
    xl <- .resolveElement(mase, x)
    yl <- .resolveElement(mase, y)
    if (!is.null(coordinate_system)) {
        g <- .maseCTgraph(mase)
        xl$layer <- .alignLayer(xl$layer, g, xl$key, coordinate_system,
                                x_col = x_col, y_col = y_col)
        yl$layer <- .alignLayer(yl$layer, g, yl$key, coordinate_system,
                                x_col = x_col, y_col = y_col)
    }
    coords <- if (xl$element_type == "points") c(x_col, y_col) else NULL
    spatialMatch(xl$layer, yl$layer, coords = coords, geom = "geometry",
                 join = join)
}
