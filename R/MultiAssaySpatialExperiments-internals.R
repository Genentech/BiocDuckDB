### =========================================================================
### MultiAssaySpatialExperiments-internals.R — MASE spatial helpers and I/O
### -------------------------------------------------------------------------

replaceSlots <- BiocGenerics:::replaceSlots

.requireDuckDBSpatial <- function(feature) {
    if (!requireNamespace("DuckDBSpatial", quietly = TRUE)) {
        stop("DuckDBSpatial required for ", feature,
             "; install with BiocManager::install('DuckDBSpatial')")
    }
}

## Capture MASE in-memory methods before BiocDuckDB overrides them.
.captureMASEMethod <- function(generic, signature = "MultiAssaySpatialExperiment") {
    if (methods::existsMethod(generic, signature))
        methods::getMethod(generic, signature)
    else
        NULL
}

.maseHasLazyLayers <- function(mase) {
    pts <- MultiAssaySpatialExperiment::spatialPoints(mase)
    shps <- MultiAssaySpatialExperiment::spatialShapes(mase)
    spmap <- MultiAssaySpatialExperiment::spatialMap(mase)
    any(vapply(as.list(pts), .isLazySpatialLayer, logical(1L))) ||
        any(vapply(as.list(shps), .isLazySpatialLayer, logical(1L))) ||
        (!is.null(spmap) && .isLazySpatialLayer(spmap))
}

.isLazySpatialLayer <- function(x) {
    inherits(x, "DuckDBDataFrame")
}

.filterSpatialLayerByInstances <- function(layer, instance_ids) {
    if (is.null(layer) || !nrow(layer) || !length(instance_ids))
        return(layer[0L, , drop = FALSE])
    if (!"instance_id" %in% colnames(layer))
        return(layer)
    ids <- as.character(instance_ids)
    if (.isLazySpatialLayer(layer)) {
        idx <- as.character(as.vector(layer[["instance_id"]])) %in% ids
        layer[idx, , drop = FALSE]
    } else {
        layer[as.character(layer[["instance_id"]]) %in% ids, , drop = FALSE]
    }
}

.filterSpatialLayersByInstances <- function(layers, instance_ids, region = NULL) {
    if (length(layers) == 0L)
        return(layers)
    out <- as.list(layers)
    nms <- if (is.null(region)) names(out) else intersect(names(out), region)
    for (nm in nms) {
        out[[nm]] <- .filterSpatialLayerByInstances(out[[nm]], instance_ids)
    }
    out
}

.propagateSpatialSubset <- function(mase, retained_by_region) {
    spmap <- MultiAssaySpatialExperiment::spatialMap(mase)
    if (is.null(spmap) || nrow(spmap) == 0L)
        return(mase)
    keep_sp <- vapply(seq_len(nrow(spmap)), function(i) {
        reg <- spmap[i, "region"]
        inst <- spmap[i, "instance_id"]
        kept <- retained_by_region[[reg]]
        !is.null(kept) && inst %in% kept
    }, logical(1L))
    spmap <- spmap[keep_sp, , drop = FALSE]
    if (nrow(spmap) == 0L)
        return(mase)
    j_list <- split(spmap[["colname"]], spmap[["assay"]])
    exps <- MultiAssayExperiment::experiments(mase)
    assay_nms <- names(exps)
    j_list_full <- lapply(assay_nms, function(a) {
        if (a %in% names(j_list)) j_list[[a]] else colnames(exps[[a]])
    })
    names(j_list_full) <- assay_nms
    MultiAssayExperiment::subsetByColumn(mase, j_list_full)
}

.spatialExtent <- function(mase) {
    xmin <- Inf; xmax <- -Inf; ymin <- Inf; ymax <- -Inf
    for (el in c(as.list(MultiAssaySpatialExperiment::spatialPoints(mase)),
                 as.list(MultiAssaySpatialExperiment::spatialShapes(mase)))) {
        if (is.null(el) || !nrow(el))
            next
        if ("x" %in% colnames(el) && "y" %in% colnames(el)) {
            xmin <- min(xmin, min(as.vector(el[["x"]]), na.rm = TRUE))
            xmax <- max(xmax, max(as.vector(el[["x"]]), na.rm = TRUE))
            ymin <- min(ymin, min(as.vector(el[["y"]]), na.rm = TRUE))
            ymax <- max(ymax, max(as.vector(el[["y"]]), na.rm = TRUE))
        }
    }
    c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Raster path references and COO labels
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom jsonlite read_json write_json
.readParquetRasterRef <- function(path, resource, ...) {
    meta_path <- file.path(path, "metadata.json")
    if (file.exists(meta_path)) {
        meta <- read_json(meta_path, simplifyVector = TRUE)
        if (!is.null(meta[["zarr_path"]]) && requireNamespace("ZarrArray", quietly = TRUE)) {
            ds <- meta[["datasets"]]
            if (is.null(ds))
                ds <- list(meta[["zarr_path"]])
            arrays <- lapply(ds, function(p) {
                ZarrArray::ZarrArray(file.path(path, p))
            })
            return(list(data = arrays, metadata = meta))
        }
        if (!is.null(meta[["file_path"]]) && requireNamespace("SpatialExperiment", quietly = TRUE)) {
            fp <- meta[["file_path"]]
            if (!file.exists(fp))
                fp <- file.path(path, fp)
            return(SpatialExperiment::SpatialImage(fp))
        }
        return(meta)
    }
    list(path = path, resource = resource)
}

.writeRasterRef <- function(x, path, name, type = c("image", "label")) {
    type <- match.arg(type)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    meta <- list(name = name, type = type)
    if (requireNamespace("SpatialExperiment", quietly = TRUE) &&
            is(x, "StoredSpatialImage")) {
        src <- SpatialExperiment::imgSource(x)
        if (length(src) && file.exists(src[[1L]])) {
            dest <- file.path(path, basename(src[[1L]]))
            if (!file.exists(dest))
                file.copy(src[[1L]], dest)
            meta[["file_path"]] <- basename(dest)
        }
    } else if (is.list(x) && length(x) && requireNamespace("ZarrArray", quietly = TRUE) &&
               all(vapply(x, function(a) inherits(a, "ZarrArray"), logical(1L)))) {
        meta[["zarr_path"]] <- path
        meta[["datasets"]] <- vapply(seq_along(x), function(i) {
            zp <- file.path(path, paste0("scale_", i))
            if (!dir.exists(zp)) dir.create(zp, recursive = TRUE)
            zp
        }, character(1L))
    } else if (is.matrix(x) || is.array(x)) {
        meta[["file_path"]] <- NA_character_
        meta[["dim"]] <- dim(x)
    }
    write_json(meta, file.path(path, "metadata.json"),
               auto_unbox = TRUE, pretty = TRUE)
    invisible(meta)
}

.matrixToCoordLabel <- function(mat, path) {
    if (!dir.exists(path))
        dir.create(path, recursive = TRUE)
    nr <- nrow(mat)
    nc <- ncol(mat)
    rows <- rep(seq_len(nr), times = nc)
    cols <- rep(seq_len(nc), each = nr)
    vals <- as.vector(mat)
    keep <- !is.na(vals) & vals != 0L
    df <- data.frame(row = rows[keep], col = cols[keep], value = vals[keep])
    DuckDBArray::writeCoordArray(df, path,
        row = list(row = seq_len(max(rows[keep], 0L))),
        col = list(col = seq_len(max(cols[keep], 0L))))
}
