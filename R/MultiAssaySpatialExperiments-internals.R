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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Raster path references and COO labels
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @importFrom jsonlite read_json write_json
.readParquetRasterRef <- function(path, resource, ...) {
    meta_path <- file.path(path, "metadata.json")
    if (!file.exists(meta_path))
        return(list(path = path, resource = resource))
    meta <- read_json(meta_path, simplifyVector = TRUE)
    nm <- meta[["name"]] %||% basename(path)
    resolve <- function(p)
        if (file.exists(p)) p else file.path(path, p)
    zp <- meta[["zarr_path"]]
    if (!is.null(zp) && !is.na(zp)) {
        zp <- resolve(zp)
        if (!requireNamespace("ZarrArray", quietly = TRUE) || !file.exists(zp)) {
            warning("raster '", nm, "' is left out: its Zarr store ", zp,
                    if (file.exists(zp)) " needs the ZarrArray package"
                    else " does not exist", call. = FALSE)
            return(NULL)
        }
        return(ZarrArray::ZarrArray(zp))
    }
    fp <- meta[["file_path"]]
    if (!is.null(fp) && !is.na(fp)) {
        fp <- resolve(fp)
        if (!requireNamespace("SpatialExperiment", quietly = TRUE) ||
                !file.exists(fp)) {
            warning("raster '", nm, "' is left out: its image file ", fp,
                    " does not exist", call. = FALSE)
            return(NULL)
        }
        return(SpatialExperiment::SpatialImage(fp))
    }
    warning("raster '", nm, "' is left out: it records no file or Zarr store",
            call. = FALSE)
    NULL
}

.writeRasterRef <- function(x, path, name, type = c("image", "label")) {
    type <- match.arg(type)
    meta <- list(name = name, type = type)
    if (requireNamespace("SpatialExperiment", quietly = TRUE) &&
            is(x, "StoredSpatialImage")) {
        src <- SpatialExperiment::imgSource(x)
        if (!length(src) || !file.exists(src[[1L]]))
            stop(type, " '", name, "' refers to an image file that does not ",
                 "exist: ", if (length(src)) src[[1L]] else "<none>",
                 call. = FALSE)
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        dest <- file.path(path, basename(src[[1L]]))
        if (!file.exists(dest))
            file.copy(src[[1L]], dest)
        meta[["file_path"]] <- basename(dest)
    } else if (is(x, "ZarrArray")) {
        zp <- DelayedArray::path(x)
        if (!length(zp) || !nzchar(zp[[1L]]))
            stop(type, " '", name, "' is a ZarrArray with no store path",
                 call. = FALSE)
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        meta[["zarr_path"]] <- normalizePath(zp[[1L]], mustWork = FALSE)
        meta[["dim"]] <- dim(x)
    } else {
        stop(type, " '", name, "' is ", paste(class(x)[1L], "in memory"),
             ", which has no file or store to refer to. Save it as an image ",
             "file (a StoredSpatialImage) or as a Zarr store (a ZarrArray), ",
             "or use a matrix for a label, which is stored as coordinates.",
             call. = FALSE)
    }
    write_json(meta, file.path(path, "metadata.json"),
               auto_unbox = TRUE, pretty = TRUE)
    invisible(meta)
}
