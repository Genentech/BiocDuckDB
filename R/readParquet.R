#' Read Parquet Representation of a Bioconductor Object
#'
#' @description
#' Reads a parquet representation of various Bioconductor objects from disk
#' using the Frictionless Data Package metadata model, reconstructing them as
#' DuckDB-backed objects. This function is the counterpart to
#' \code{\link{writeParquet}} and supports reading \code{SummarizedExperiment},
#' \code{RangedSummarizedExperiment}, \code{SingleCellExperiment}, and
#' \code{MultiAssayExperiment} objects that have been written to parquet format
#' with comprehensive metadata.
#'
#' @param path A character string specifying the path to the directory
#' containing the parquet files. The directory must contain a
#' \code{datapackage.json} file written by \code{\link{writeParquet}}.
#' @param package A list containing the parsed \code{datapackage.json}
#' contents. Defaults to reading \code{file.path(path, "datapackage.json")}.
#' Contains the package-level \code{model} field and a \code{resources} list
#' describing each parquet dataset.
#' @param model A character string identifying the package schema — which
#' container class to reconstruct. Defaults to \code{package[["model"]]}.
#' Supported values: \code{"summarized_experiment"},
#' \code{"ranged_summarized_experiment"}, \code{"single_cell_experiment"},
#' \code{"experiment_list"}, \code{"multi_assay_experiment"},
#' \code{"multi_assay_spatial_experiment"}. When \code{NULL} (absent from the
#' package JSON), all resources are returned as a \code{SimpleList} with no
#' imposed schema.
#' @param ... Additional arguments passed to internal helper functions.
#'
#' @return
#' A Bioconductor object of the appropriate class, reconstructed using DuckDB
#' backend classes:
#' \itemize{
#'   \item \code{\linkS4class{DuckDBMatrix}} objects for assay data
#'   \item \code{\linkS4class{DuckDBDataFrame}} objects for row and column metadata
#'   \item \code{\linkS4class{DuckDBGRanges}} objects for genomic ranges
#'   \item \code{\linkS4class{DuckDBGRangesList}} objects for grouped genomic ranges
#'   \item \code{\linkS4class{DuckDBSelfHits}} objects for graph edge lists
#'   \item \code{\linkS4class{DuckDBDualSubset}} objects for pairwise graphs (wrapping \code{DuckDBSelfHits})
#' }
#'
#' @details
#' This function reads parquet files that were created by
#' \code{\link{writeParquet}} using the Frictionless Data Package metadata model
#' and reconstructs the original Bioconductor object structure. This function
#' supports different object types through a switch statement that routes to
#' appropriate helper functions:
#'
#' \strong{\code{SummarizedExperiment} objects:} Reads feature data, sample data,
#' and assay data from separate parquet files, reconstructing the object with
#' DuckDB-backed components. For ranged experiments, genomic coordinates are
#' reconstructed as \code{DuckDBGRanges} objects.
#'
#' \strong{\code{SingleCellExperiment} objects:} Extends \code{SummarizedExperiment}
#' functionality by also reading reduced dimensions, loadings, alternative
#' experiments, row / column tables, and row / column pairings from their
#' respective parquet files.
#'
#' \strong{\code{MultiAssayExperiment} objects:} Reads multiple experiments
#' along with sample data and sample mapping information, reconstructing the
#' multi-experiment structure.
#'
#' This function uses the Frictionless Data Package metadata to understand the
#' structure of the data, including resource schemas with field definitions,
#' primary keys, and semantic meanings (sequence_name, start, end, strand for
#' genomic data) encoded in the \code{datapackage.json} file.
#'
#' @section Supported Object Types:
#' \describe{
#'   \item{\code{GenomicRanges}}{
#'     Read with genomic coordinates (\code{seqnames}, \code{start},
#'     \code{end}, \code{strand}) and optional metadata columns.
#'   }
#'   \item{\code{GenomicRangesList}}{
#'     Read with genomic coordinates and metadata columns.
#'   }
#'   \item{\code{DuckDBSelfHits}}{
#'     Graph edge lists with \code{from}, \code{to} columns and optional
#'     metadata. Node count (\code{nnode}) is stored in the schema
#'     \code{graphEdges} property.
#'   }
#'   \item{\code{SummarizedExperiment}}{
#'     Feature metadata from \code{features/}, sample metadata from
#'     \code{samples/}, and assays from flat \code{assay_<name>/} directories.
#'     Any complex objects stored in \code{metadata()} are written to
#'     \code{unbound_<name>/} directories and restored on read.
#'   }
#'   \item{\code{RangedSummarizedExperiment}}{
#'     As \code{SummarizedExperiment}, with \code{rowRanges} reconstructed as a
#'     \code{DuckDBGRanges} from \code{features/}.
#'   }
#'   \item{\code{SingleCellExperiment}}{
#'     Extends \code{SummarizedExperiment} with: reduced dimensions from
#'     \code{sample_embeddings/}; row loadings from \code{feature_embeddings/};
#'     alternative experiments from \code{experiment_<name>/} directories;
#'     row/column tables from flat \code{feature_table_<name>/} and
#'     \code{sample_table_<name>/} directories; row/column pairwise graphs from
#'     flat \code{feature_graph_<name>/} and \code{sample_graph_<name>/}
#'     directories.
#'   }
#'   \item{\code{MultiAssayExperiment}}{
#'     Experiments written directly to root as \code{experiment_<name>/}
#'     directories (flattened). Each \code{SummarizedExperiment} has its own
#'     \code{datapackage.json} (\code{layout = "nested_experiment"}).
#'     Also includes \code{subjects} (subject metadata) and \code{sample_map}
#'     (subject-to-sample mapping) as \code{unbound} resources.
#'   }
#'   \item{\code{MultiAssaySpatialExperiment}}{
#'     Extends \code{MultiAssayExperiment} with spatial points from flat
#'     \code{sample_points_<name>/} directories, shapes from flat
#'     \code{sample_shapes_<name>/} directories, image metadata from
#'     \code{img_data/}, and spatial mapping from \code{spatial_map/}.
#'     Requires the \code{MultiAssaySpatialExperiment} package.
#'   }
#' }
#'
#' @section Frictionless Data Package Metadata:
#' This function reads data packages that follow the Frictionless Data Package
#' specification (extended by the BiocDuckDB profile), parsing the
#' \code{datapackage.json} metadata file. Dispatch operates at two levels:
#'
#' \strong{Package level} (root \code{datapackage.json}):
#' \itemize{
#'   \item \code{model} — which container class to reconstruct; absent means
#'     \code{SimpleList} fallback.
#'   \item \code{annotations} — non-relational metadata from \code{metadata(x)}
#'     (scalars, vectors, 1-D arrays, JSON-safe nested lists). Tabular items
#'     are referenced via \code{parquet_ref} stubs pointing at \code{unbound}
#'     Parquet sidecars; mixed nested lists use a \code{nested_mapping} wrapper.
#' }
#'
#' \strong{Resource level} (entries in \code{resources} array):
#' \itemize{
#'   \item \code{dimension} — biological axis (\code{"feature"},
#'     \code{"sample"}, \code{"crossed"}, \code{"unbound"}).
#'   \item \code{layout} — physical storage pattern (e.g.,
#'     \code{"data_frame"}, \code{"coord_array"}, \code{"graph_edges"},
#'     \code{"nested_experiment"}); routes each resource to the
#'     appropriate low-level reader via \code{.readParquetResource}.
#'   \item \code{schema} — field definitions, primary keys, and BiocDuckDB
#'     annotations (\code{genomicCoords}, \code{graphEdges}, \code{arrayItem}).
#' }
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{writeParquet}} for writing Bioconductor objects to parquet
#'   \item \code{\linkS4class{DuckDBMatrix}} for matrix storage
#'   \item \code{\linkS4class{DuckDBDataFrame}} for metadata storage
#'   \item \code{\linkS4class{DuckDBGRanges}} for genomic ranges
#'   \item \code{\linkS4class{DuckDBGRangesList}} for grouped genomic ranges
#'   \item \code{\linkS4class{DuckDBSelfHits}} for graph edge lists
#' }
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(rpois(500, 5), 25, 20)),
#'     colData = S4Vectors::DataFrame(sample = rep(c("A", "B"), each = 10)))
#' tmpdir <- tempfile()
#' writeParquet(se, tmpdir)
#' se2 <- readParquet(tmpdir)
#' class(se2)
#' dim(se2)
#' unlink(tmpdir, recursive = TRUE)
#'
#' @include fieldtypes.R
#'
#' @export
#' @rdname readParquet
readParquet <-
function(path,
         package = .readDatapackage(path),
         model = package[["model"]],
         ...)
{
    if (is.null(model))
        return(.readParquetSimpleList(path, package, ...))
    switch(model,
           "summarized_experiment"          =,
           "ranged_summarized_experiment"   =,
           "single_cell_experiment"         = .readParquetSE(path, package, ...),
           "experiment_list"                = .readParquetExps(path, package, ...),
           "multi_assay_experiment"         = .readParquetMAE(path, package, ...),
           "multi_assay_spatial_experiment" = .readParquetMASE(path, package, ...),
           stop("unsupported model: ", model))
}

#' @importFrom jsonlite fromJSON read_json
#' @importFrom DBI dbGetQuery
#' @importFrom DuckDBDataFrame acquireDuckDBConn configureCloud
.readDatapackage <- function(path) {
    dp <- file.path(path, "datapackage.json")
    remote <- is.character(path) && length(path) == 1L && !is.na(path) &&
        grepl("^(s3|gs|gcs|az|azure|abfss|r2|http|https)://", path)
    if (remote) {
        conn <- acquireDuckDBConn()
        configureCloud(conn)
        res <- dbGetQuery(conn, sprintf(
            "SELECT content FROM read_text('%s')", gsub("'", "''", dp)))
        if (nrow(res) == 0L || is.na(res[["content"]][[1L]])) {
            stop("could not read remote datapackage.json at '", dp, "'")
        }
        return(fromJSON(res[["content"]][[1L]], simplifyVector = TRUE,
                        simplifyDataFrame = FALSE, simplifyMatrix = FALSE))
    }
    read_json(dp, simplifyVector = TRUE, simplifyDataFrame = FALSE,
              simplifyMatrix = FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

.schema_keycols <- function(schema) {
    keycol <- schema[["sortOrder"]][[1L]][["field"]] %||% schema[["primaryKey"]]
    if (is.null(keycol)) {
        keycol <- unlist(lapply(schema[["foreignKeys"]], `[[`, "fields"))
    }
    keycol
}

.schema_datacols <- function(schema) {
    all_fields <- sapply(schema[["fields"]], `[[`, "name")
    keycol <- .schema_keycols(schema)
    pkey <- schema[["primaryKey"]]
    exclude <- c(keycol,
                 if (!is.null(pkey) && !identical(pkey, keycol)) pkey,
                 .schema_partitions(schema),
                 unlist(schema[["genomicCoords"]]),
                 unlist(schema[["graphEdges"]]))
    setdiff(all_fields, exclude)
}

## Factor levels recorded for the given data columns, in DuckDBTable's
## 'collevels' shape (named list of list(levels, ordered)). Only the flat-table
## factor shape (categories = character vector of levels, from
## fieldtypes.R:.fieldtype.factor) is restored here; the array/keycol shape
## (categories = list of {value, label} structs) is left to the array reader.
.schema_collevels <- function(schema, datacols) {
    fields <- schema[["fields"]]
    if (length(fields) == 0L) {
        return(NULL)
    }
    names(fields) <- sapply(fields, `[[`, "name")
    out <- list()
    for (nm in intersect(datacols, names(fields))) {
        cats <- fields[[nm]][["categories"]]
        if (is.null(cats) || length(cats) == 0L) {
            next
        }
        first <- cats[[1L]]
        if (is.list(first) && !is.null(first[["value"]])) {
            next  # array/keycol {value,label} shape, not a flat factor
        }
        out[[nm]] <- list(levels = as.character(unlist(cats)),
                          ordered = isTRUE(fields[[nm]][["categoriesOrdered"]]))
    }
    if (length(out) == 0L) NULL else out
}

.schema_genomic <- function(schema, role) {
    schema[["genomicCoords"]][[role]]
}

.schema_graph <- function(schema, role) {
    schema[["graphEdges"]][[role]]
}

.schema_partitions <- function(schema) {
    schema[["partitioning"]]
}

.isMetadataStub <- function(x) {
    is.list(x) && !is.null(x[["__type__"]]) && is.character(x[["__type__"]])
}

.referencedMetadataResources <- function(annotations) {
    .collectRefs <- function(x) {
        if (!.isMetadataStub(x)) {
            if (is.list(x) && !is.data.frame(x))
                return(unlist(lapply(x, .collectRefs), use.names = FALSE))
            return(character(0L))
        }
        if (identical(x[["__type__"]], "parquet_ref"))
            return(x[["resource"]])
        if (identical(x[["__type__"]], "nested_mapping")) {
            nms <- setdiff(names(x), "__type__")
            return(unlist(lapply(x[nms], .collectRefs), use.names = FALSE))
        }
        character(0L)
    }
    unique(unlist(lapply(annotations, .collectRefs), use.names = FALSE))
}

.deserializeMetadataValue <- function(value, resources_by_name, path, ...) {
    if (!.isMetadataStub(value)) {
        if (is.list(value) && !is.data.frame(value)) {
            if (any(vapply(value, .isMetadataStub, logical(1L)))) {
                return(lapply(value, .deserializeMetadataValue,
                              resources_by_name = resources_by_name,
                              path = path, ...))
            }
        }
        return(value)
    }
    switch(value[["__type__"]],
           parquet_ref = {
               res <- resources_by_name[[value[["resource"]]]]
               if (is.null(res)) {
                   stop("Metadata resource '", value[["resource"]],
                        "' not found in datapackage.json")
               }
               .readParquetResource(path, res, ...)
           },
           nested_mapping = {
               children <- value[setdiff(names(value), "__type__")]
               lapply(children, .deserializeMetadataValue,
                      resources_by_name = resources_by_name,
                      path = path, ...)
           },
           value)
}

.deserializeMetadata <- function(annotations, resources, path, ...) {
    annotations <- annotations %||% list()
    resources_by_name <- setNames(resources,
                                  vapply(resources, `[[`, character(1L), "name"))
    metadata <- lapply(annotations, .deserializeMetadataValue,
                       resources_by_name = resources_by_name,
                       path = path, ...)
    referenced <- .referencedMetadataResources(annotations)
    for (res in .filterResources(resources, "unbound")) {
        nm <- res[["name"]]
        path_nm <- res[["path"]]
        if (is.null(metadata[[nm]]) && !(nm %in% referenced) &&
            startsWith(path_nm, "unbound_")) {
            metadata[[nm]] <- .readParquetResource(path, res, ...)
        }
    }
    metadata
}

.filterResources <- function(resources, dimension, layout = NULL) {
    if (is.null(layout)) {
        fun <- function(x) isTRUE(x[["dimension"]] %in% dimension)
    } else {
        fun <- function(x) {
            isTRUE(x[["dimension"]] %in% dimension) &&
            isTRUE(x[["layout"]] %in% layout)
        }
    }
    Filter(fun, resources)
}

#' @importFrom S4Vectors SimpleList
.readParquetSimpleList <- function(path, package, ...) {
    resources <- package[["resources"]]
    out <- lapply(resources, function(r) .readParquetResource(path, r, ...))
    names(out) <- sapply(resources, `[[`, "name")
    SimpleList(out)
}

.readParquetResource <- function(path, resource, ...) {
    fullpath <- file.path(path, resource[["path"]])
    switch(resource[["layout"]],
           "data_frame"            =,
           "embedding_table"       =,
           "nested_data_frame"     =,
           "spatial_points"        =,
           "spatial_shapes"        = .readParquetSpatialDataFrame(fullpath, resource, ...),
           "transposed_data_frame" = .readParquetTransposedDataFrame(fullpath, resource, ...),
           "coord_array"           = .readParquetArray(fullpath, resource, ...),
           "genomic_ranges"        = .readParquetGenomicRanges(fullpath, resource, ...),
           "genomic_ranges_list"   = .readParquetGenomicRangesList(fullpath, resource, ...),
           "spatial_raster_ref"    = .readParquetRasterRef(fullpath, resource, ...),
           "spatial_label_coord"   = .readParquetArray(fullpath, resource, ...),
           "graph_edges"           = .readParquetGraphEdges(fullpath, resource, ...),
           stop("unsupported layout: ", resource[["layout"]]))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DataFrame objects
###

#' @importFrom DuckDBDataFrame DuckDBDataFrame
.readParquetDataFrame <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    if (length(keycol) > 1L) {
        keycol <- keycol[1L]
    }
    DuckDBDataFrame(path, datacols = datacols, keycol = keycol,
                    collevels = .schema_collevels(resource[["schema"]], datacols))
}

#' @importFrom DuckDBDataFrame DuckDBDataFrame acquireDuckDBConn
.readParquetSpatialDataFrame <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    if (resource[["layout"]] == "spatial_shapes" &&
        requireNamespace("DuckDBSpatial", quietly = TRUE)) {
        DuckDBSpatial::enableGeoParquetConversion(acquireDuckDBConn())
    }
    if (length(keycol) > 1L)
        keycol <- keycol[1L]
    DuckDBDataFrame(path, datacols = datacols, keycol = keycol,
                    collevels = .schema_collevels(resource[["schema"]], datacols))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TransposedDataFrame objects
###

#' @importFrom DuckDBDataFrame DuckDBDataFrame
.readParquetTransposedDataFrame <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    if (length(keycol) > 1L) {
        keycol <- keycol[1L]
    }
    t(DuckDBDataFrame(path, datacols = datacols, keycol = keycol))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRanges objects
###

#' @importFrom DuckDBGRanges DuckDBGRanges
.readParquetGenomicRanges <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    schema <- resource[["schema"]]
    DuckDBGRanges(path,
                  seqnames = .schema_genomic(schema, "seqname"),
                  start    = .schema_genomic(schema, "start"),
                  end      = .schema_genomic(schema, "end"),
                  strand   = .schema_genomic(schema, "strand"),
                  mcols    = datacols,
                  keycol   = keycol)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRangesList objects
###

#' @importFrom DuckDBGRanges DuckDBGRangesList
#' @importFrom S4Vectors mcols mcols<-
.readParquetGenomicRangesList <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    schema <- resource[["schema"]]
    DuckDBGRangesList(path,
                      seqnames = .schema_genomic(schema, "seqname"),
                      start    = .schema_genomic(schema, "start"),
                      end      = .schema_genomic(schema, "end"),
                      strand   = .schema_genomic(schema, "strand"),
                      mcols    = datacols,
                      keycol   = keycol)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Graph edges (SelfHits) objects
###

#' @importFrom DuckDBDataFrame DuckDBSelfHits
.readParquetGraphEdges <-
function(path,
         resource,
         datacols = .schema_datacols(resource[["schema"]]),
         keycol = .schema_keycols(resource[["schema"]]),
         ...)
{
    schema <- resource[["schema"]]
    DuckDBSelfHits(path,
                   from = .schema_graph(schema, "from"),
                   to = .schema_graph(schema, "to"),
                   nnode = .schema_graph(schema, "nnode"),
                   mcols = datacols,
                   keycol = keycol)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Array objects
###

#' @importFrom DuckDBArray DuckDBArray DuckDBMatrix
#' @importFrom stats setNames
.readParquetArray <-
function(path,
         resource,
         keycols = NULL,
         dimtbls = NULL, ...)
{
    if (is.null(keycols)) {
        keycols <- .schema_keycols(resource[["schema"]])
        fields <- resource[["schema"]][["fields"]]
        names(fields) <- sapply(fields, `[[`, "name")
        keycols <- lapply(fields[keycols], function(x) {
            if (is.null(x[["categories"]])) {
                NULL
            } else {
                sapply(x[["categories"]], function(y) {
                    sort(setNames(y[["value"]], y[["label"]]))
                })
            }
        })
    }
    # For arrays, the datacol is the single value column (all datacols from schema)
    datacols <- .schema_datacols(resource[["schema"]])
    datacol <- datacols[1L]  # Arrays have a single data column
    
    if (length(keycols) == 2L) {
        DuckDBMatrix(path, datacol = datacol, keycols = keycols, dimtbls = dimtbls)
    } else {
        DuckDBArray(path, datacol = datacol, keycols = keycols, dimtbls = dimtbls)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SummarizedExperiment objects
###

#' @importFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom IRanges DataFrameList
#' @importFrom S4Vectors endoapply
#' @importFrom SingleCellExperiment SingleCellExperiment colPair<- rowPair<-
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment SummarizedExperiment
.readParquetSE <- function(path, package, ...) {
    resources <- package[["resources"]]
    names(resources) <- sapply(package[["resources"]], `[[`, "name")

    # Dimension Resources
    feature_res <- resources[["features"]]
    sample_res <- resources[["samples"]]

    # Dimension Tables
    fun <- function(res) {
        schema <- res[["schema"]]
        keycol <- .schema_keycols(schema)
        pkey <- schema[["primaryKey"]]
        name <- if (!identical(pkey, keycol)) pkey else NULL
        partitions <- .schema_partitions(schema)
        df <- DuckDBDataFrame(file.path(path, res[["path"]]),
                              datacols = c(name, partitions),
                              keycol = keycol)
        df <- as.data.frame(df, optional = TRUE)
        if (length(name) > 0L) {
            rownames(df) <- df[[1L]]
            df <- df[-1L]
        }
        setNames(list(df), keycol)
    }
    dimtbls <- new.env(parent = emptyenv())
    dimtbls[["dimtbls"]] <- as(c(fun(feature_res), fun(sample_res)), "DataFrameList")
    dimkeycols <- lapply(dimtbls[["dimtbls"]],
                         function(x) setNames(seq_len(nrow(x)), rownames(x)))
    dimtbls[["dimtbls"]] <- endoapply(dimtbls[["dimtbls"]], function(x) {
        rownames(x) <- NULL
        x
    })

    # Feature Data
    keycol <- dimkeycols[1L]
    feature_data <- NULL
    feature_ranges <- NULL
    if (isTRUE(feature_res[["layout"]] %in% c("genomic_ranges", "genomic_ranges_list"))) {
        feature_ranges <- .readParquetResource(path, feature_res, keycol = keycol)
    } else {
        feature_data <- .readParquetResource(path, feature_res, keycol = keycol)
    }

    # Sample Data
    keycol <- dimkeycols[2L]
    samples <- .readParquetResource(path, sample_res, keycol = keycol)

    # Assays
    assay_res <- .filterResources(resources, "crossed",
                                  c("coord_array", "data_frame",
                                    "transposed_data_frame"))
    assays <- lapply(assay_res, function(res) {
        # Assays are written transposed, so we need to transpose them back
        t(.readParquetResource(path, res, keycol = rev(dimkeycols),
                               dimtbls = dimtbls))
    })
    names(assays) <- sapply(assay_res, `[[`, "name")

    # Metadata
    metadata <- .deserializeMetadata(package[["annotations"]], resources, path, ...)

    # SummarizedExperiment
    if (isTRUE(package[["model"]] %in% c("ranged_summarized_experiment", "summarized_experiment"))) {
        se <- SummarizedExperiment(assays,
                                   rowData = feature_data,
                                   rowRanges = feature_ranges,
                                   colData = samples,
                                   metadata = metadata)
    } else if (isTRUE(package[["model"]] == "single_cell_experiment")) {
        # Reduced Dimensions
        rdim_res <- .filterResources(resources, "sample", "embedding_table")
        rdims <- if (length(rdim_res)) {
            as.list(.readParquetResource(path, rdim_res[[1L]], keycol = dimkeycols[2L]))
        } else {
            list()
        }

        # Alternative Experiments
        altexp_res <- .filterResources(resources, "crossed",
                                       "nested_experiment")
        alts <- if (length(altexp_res)) {
            as.list(.readParquetExps(path, altexp_res, ...))
        } else {
            list()
        }

        # SingleCellExperiment
        se <- SingleCellExperiment(assays,
                                   rowData = feature_data,
                                   rowRanges = feature_ranges,
                                   colData = samples,
                                   reducedDims = rdims,
                                   altExps = alts,
                                   mainExpName = package[["main_exp_name"]],
                                   metadata = metadata)

        # Row Loadings
        load_res <- .filterResources(resources, "feature", "embedding_table")
        if (length(load_res)) {
            rowLoadings(se) <- as.list(.readParquetResource(path, load_res[[1L]],
                                                            keycol = dimkeycols[1L]))
        }

        # Row Tables
        for (res in .filterResources(resources, "feature", "nested_data_frame")) {
            rowTable(se, res[["name"]]) <- .readParquetResource(path, res,
                                                                keycol = dimkeycols[1L])
        }

        # Column Tables
        for (res in .filterResources(resources, "sample", "nested_data_frame")) {
            colTable(se, res[["name"]]) <- .readParquetResource(path, res,
                                                                keycol = dimkeycols[2L])
        }

        # Row Pairs
        for (res in .filterResources(resources, "feature", "graph_edges")) {
            rowPair(se, res[["name"]]) <- .readParquetResource(path, res)
        }

        # Column Pairs
        for (res in .filterResources(resources, "sample", "graph_edges")) {
            colPair(se, res[["name"]]) <- .readParquetResource(path, res)
        }
    }

    se
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ExperimentList objects
###

#' @importFrom MultiAssayExperiment ExperimentList
.readParquetExps <- function(path, resources, ...) {
    exps <- lapply(resources, function(res) {
        if (isTRUE(res[["layout"]] == "nested_experiment")) {
            readParquet(file.path(path, res[["path"]]), ...)
        } else { # array-like object
            .readParquetResource(path, res, ...)
        }
    })
    names(exps) <- sapply(resources, `[[`, "name")
    ExperimentList(exps)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MultiAssayExperiment objects
###

#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom stats setNames
.readParquetMAE <- function(path, package, ...) {
    resources <- package[["resources"]]
    names(resources) <- sapply(package[["resources"]], `[[`, "name")

    # Experiments
    exps_res <- .filterResources(resources, dimension = "crossed")
    experiments <- .readParquetExps(path, exps_res, ...)

    # Subject Data
    subject_res <- resources[["subjects"]]
    subject_schema <- subject_res[["schema"]]
    index <- .schema_keycols(subject_schema)
    pkey <- subject_schema[["primaryKey"]]
    name <- if (!is.null(pkey) && !identical(pkey, index)) pkey else NULL
    fullpath <- file.path(path, subject_res[["path"]])
    keycol <- .readParquetDataFrame(fullpath, resource = subject_res,
                                    datacols = name, keycol = index)
    keycol <- as.list(as.data.frame(keycol, optional = TRUE))
    subjects <- .readParquetDataFrame(fullpath, resource = subject_res,
                                      keycol = keycol)

    # Sample Map
    sample_map_res <- resources[["sample_map"]]
    index <- .schema_keycols(sample_map_res[["schema"]])
    fullpath <- file.path(path, sample_map_res[["path"]])
    sample_map <- .readParquetDataFrame(fullpath, resource = sample_map_res,
                                        keycol = index)
    sample_map <- as.data.frame(sample_map, optional = TRUE)
    sample_map[[1L]] <- factor(sample_map[[1L]], levels = names(experiments))

    # Metadata
    metadata <- .deserializeMetadata(package[["annotations"]], resources, path, ...)

    # MultiAssayExperiment
    MultiAssayExperiment(experiments,
                         colData = subjects,
                         sampleMap = sample_map,
                         metadata = metadata)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MultiAssaySpatialExperiment objects
### -------------------------------------------------------------------------

#' @importFrom jsonlite read_json
#' @importFrom S4Vectors DataFrame
.readParquetMASE <- function(path, package, ...) {
    if (!requireNamespace("MultiAssaySpatialExperiment", quietly = TRUE)) {
        stop("package 'MultiAssaySpatialExperiment' is required to read ",
             "multi_assay_spatial_experiment; install it with ",
             "BiocManager::install(\"MultiAssaySpatialExperiment\")")
    }
    resources <- package[["resources"]]
    names(resources) <- sapply(package[["resources"]], `[[`, "name")

    # Base MAE object
    mae <- .readParquetMAE(path, package, ...)

    # Points (optional) - keep lazy DuckDBDataFrame
    points_res <- .filterResources(resources, "sample", "spatial_points")
    if (length(points_res)) {
        points_list <- lapply(points_res, function(r) {
            .readParquetResource(path, r, ...)
        })
        names(points_list) <- sub("^sample_points_", "",
            sapply(points_res, `[[`, "name"))
        points <- MultiAssaySpatialExperiment::PointsLayerList(points_list)
    } else {
        points <- MultiAssaySpatialExperiment::PointsLayerList()
    }

    # Shapes (optional) - keep lazy DuckDBDataFrame
    shapes_res <- .filterResources(resources, "sample", "spatial_shapes")
    if (length(shapes_res)) {
        shapes_list <- lapply(shapes_res, function(r) {
            .readParquetResource(path, r, ...)
        })
        names(shapes_list) <- sub("^sample_shapes_", "",
            sapply(shapes_res, `[[`, "name"))
        shapes <- MultiAssaySpatialExperiment::ShapesLayerList(shapes_list)
    } else {
        shapes <- MultiAssaySpatialExperiment::ShapesLayerList()
    }

    # Images (optional) - raster path references
    images <- MultiAssaySpatialExperiment::RasterLayerList()
    img_res <- .filterResources(resources, "sample", "spatial_raster_ref")
    if (length(img_res)) {
        images_list <- lapply(img_res, function(r) {
            .readParquetRasterRef(file.path(path, r[["path"]]), r, ...)
        })
        names(images_list) <- sub("^sample_images_", "", sapply(img_res, `[[`, "name"))
        images <- MultiAssaySpatialExperiment::RasterLayerList(images_list)
    } else if (!is.null(resources[["images"]])) {
        fullpath <- file.path(path, resources[["images"]][["path"]])
        if (dir.exists(fullpath)) {
            subdirs <- list.dirs(fullpath, recursive = FALSE, full.names = TRUE)
            if (length(subdirs)) {
                images_list <- lapply(subdirs, function(d) {
                    .readParquetRasterRef(d, list(path = basename(d)), ...)
                })
                names(images_list) <- basename(subdirs)
                images <- MultiAssaySpatialExperiment::RasterLayerList(images_list)
            }
        }
    }

    # Labels (optional) - raster refs or COO coord arrays
    labels <- MultiAssaySpatialExperiment::RasterLayerList()
    lab_coord <- .filterResources(resources, "sample", "spatial_label_coord")
    lab_ref <- .filterResources(resources, "sample", "spatial_raster_ref")
    lab_ref <- lab_ref[grep("^sample_labels_", sapply(lab_ref, `[[`, "name"))]
    if (length(lab_coord)) {
        labels_list <- lapply(lab_coord, function(r) {
            .readParquetArray(file.path(path, r[["path"]]), r, ...)
        })
        names(labels_list) <- sub("^sample_labels_", "", sapply(lab_coord, `[[`, "name"))
        labels <- MultiAssaySpatialExperiment::RasterLayerList(labels_list)
    }
    if (length(lab_ref)) {
        ref_list <- lapply(lab_ref, function(r) {
            .readParquetRasterRef(file.path(path, r[["path"]]), r, ...)
        })
        names(ref_list) <- sub("^sample_labels_", "", sapply(lab_ref, `[[`, "name"))
        if (length(labels))
            labels[names(ref_list)] <- ref_list
        else
            labels <- MultiAssaySpatialExperiment::RasterLayerList(ref_list)
    }

    # imgData (optional) - must materialize
    img_data <- NULL
    if (!is.null(resources[["img_data"]])) {
        fullpath <- file.path(path, resources[["img_data"]][["path"]])
        img_data <- .readParquetDataFrame(fullpath,
                                          resource = resources[["img_data"]],
                                          ...)
        img_data <- as(img_data, "DFrame")

        # Reconstruct data column with SpatialImage objects from materialized files
        if ("image_file" %in% colnames(img_data)) {
            # Reconstruct StoredSpatialImage objects from saved PNG files
            img_objs <- lapply(seq_len(nrow(img_data)), function(i) {
                img_relpath <- img_data[["image_file"]][i]

                if (is.na(img_relpath)) {
                    return(NULL)
                }

                # Convert relative path to absolute path
                img_fullpath <- file.path(path, img_relpath)

                if (file.exists(img_fullpath)) {
                    # Create StoredSpatialImage pointing to the materialized PNG
                    if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
                        SpatialExperiment::SpatialImage(img_fullpath)
                    } else {
                        # Fallback: store path as character
                        img_fullpath
                    }
                } else {
                    warning("Image file not found: ", img_fullpath)
                    NULL
                }
            })

            # Remove image_file column and add reconstructed data column
            img_data[["image_file"]] <- NULL
            img_data[["data"]] <- img_objs
        } else {
            # Add data column with NULL values for compatibility if missing
            if (!"data" %in% colnames(img_data)) {
                img_data[["data"]] <- rep(list(NULL), nrow(img_data))
            }
        }
    }

    # spatialMap (optional) - materialize for MASE validation
    spatial_map <- NULL
    if (!is.null(resources[["spatial_map"]])) {
        spatial_map_res <- resources[["spatial_map"]]
        index <- .schema_keycols(spatial_map_res[["schema"]])
        fullpath <- file.path(path, spatial_map_res[["path"]])
        x <- .readParquetDataFrame(fullpath, resource = spatial_map_res,
                                   keycol = index, ...)
        spatial_map <- as.data.frame(x, optional = TRUE)
        spatial_map <- DataFrame(spatial_map, check.names = FALSE)
    }

    MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
        experiments = experiments(mae),
        colData = colData(mae),
        sampleMap = sampleMap(mae),
        images = images,
        labels = labels,
        points = points,
        shapes = shapes,
        imgData = img_data,
        spatialMap = spatial_map,
        metadata = metadata(mae)
    )
}
