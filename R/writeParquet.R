#' Write Parquet Representation of a Bioconductor Object
#'
#' @description
#' An \code{arrow::write_dataset} wrapper function to write the Parquet
#' representation of various R objects using the Frictionless
#' [Data Package](https://datapackage.org) metadata model. This function
#' supports multiple object types including arrays, data frames, genomic
#' ranges, and more, converting them into efficient Parquet format for
#' analysis and storage with comprehensive metadata descriptions.
#'
#' @param x The object to write. Supported types include:
#' \itemize{
#'   \item Array-like objects - converted to coordinate (long) format
#'   \item \code{data.frame} objects - written directly with optional row names
#'   \item \code{DataFrame} objects - written with column type information
#'   \item \code{TransposedDataFrame} objects - transposed before writing
#'   \item \code{GenomicRanges} objects - genomic coordinates with metadata
#'   \item \code{GenomicRangesList} objects - lists of genomic ranges
#'   \item \code{SelfHits} objects - graph edge lists with node count
#'   \item \code{Assays} objects - multi-assay arrays
#'   \item \code{SummarizedExperiment} objects - multi-assay summarized experiments
#'   \item \code{RangedSummarizedExperiment} objects - multi-assay genomic experiments
#'   \item \code{SingleCellExperiment} objects - single-cell genomic experiments
#'   \item \code{ExperimentList} objects - multi-experiment lists
#'   \item \code{MultiAssayExperiment} objects - multi-experiment genomic studies
#' }
#' @param path A character string specifying the path where the data will be
#' written. The path will be created if it doesn't exist.
#' @param indexcol For \code{data.frame} and \code{Vector} objects, a character
#' string specifying the column name for the integer index column in the
#' resulting table. Defaults to "__index__". Set to \code{NULL} to omit this
#' column.
#' @param keycol For \code{data.frame} and \code{Vector} objects, a character
#' string specifying the column name for the row names (data frames) / names
#' (genomic ranges) column in the resulting table. Defaults to "__name__".
#' Set to \code{NULL} to omit this column.
#' @param dimtbl For \code{data.frame} and \code{Vector} objects, an optional
#' \code{data.frame} containing a dimension table with partitioning information.
#' @param indexcols For array-like objects, a character vector of column names
#' for the integer index columns in the resulting table. Defaults to the
#' \code{names(dimnames(x))} or \code{sprintf("index\%d", seq_along(dim(x)))}
#' if dimension names are not available.
#' @param indexrefs For array-like objects, an optional list of foreign key
#' references for the index columns in the schema. Defaults to \code{NULL}.
#' @param datacol For array-like objects, a character string specifying the
#' column name containing the array values in the resulting table. Defaults to
#' "value".
#' @param grid For array-like objects, an optional \link[S4Arrays]{ArrayGrid}
#' to use for partitioning array-like objects. If provided, the array will be
#' split into chunks and written as separate files. Defaults to
#' \code{defaultAutoGrid(COO_SparseArray(dim(x)))}.
#' @param grid_suffix For array-like objects, a character string to append to
#' the partitioning directories if the grid contains more than one cell.
#' Defaults to "_group".
#' @param BPPARAM For array-like objects, a
#' \link[BiocParallel]{BiocParallelParam} object to use for parallel processing
#' when \code{grid} contains multiple cells. Defaults to
#' \code{getAutoBPPARAM()}.
#' @param name A character string specifying the name of the resource. Defaults
#' to the basename of the \code{path} argument.
#' @param class A character string specifying the class of the resource.
#' @param package A list containing the name of the data package and a list of
#' resources. Only used by complex objects such as \code{SummarizedExperiment}
#' and \code{MultiAssayExperiment}.
#' @param ... Additional arguments to pass to \code{arrow::write_dataset},
#' such as \code{format} (e.g., "parquet", "csv"), \code{compression}, etc.
#'
#' @return Invisibly returns the \code{path} argument, allowing for function
#' chaining.
#'
#' @details
#' This function provides specialized handling for different object types:
#'
#' \strong{Array-like objects:} Converts multi-dimensional arrays into a
#' coordinate (long) format where each non-zero element is represented as a
#' row with columns for each dimension and the value. For sparse arrays, only
#' non-zero elements are written, making it efficient for sparse data. When a
#' grid is provided with multiple cells, the array is partitioned and each
#' partition is written to a separate subdirectory.
#'
#' \strong{\code{data.frame} objects}: Writes data frames directly with optional
#' row names and dimension lookup tables for partitioning information.
#'
#' \strong{\code{DataFrame} objects:} Writes data frames with schema properties
#' that describe column types and semantics. Column metadata (\code{mcols}) is
#' not preserved; use object-level \code{metadata()} for annotations.
#'
#' \strong{\code{TransposedDataFrame} objects:} Writes transposed data frames
#' with optional row names and dimension lookup tables for partitioning
#' information.
#'
#' \strong{\code{GenomicRanges} objects:} Writes genomic coordinates with proper
#' column naming and schema properties (genomicCoords) that describe the mapping
#' of genomic coordinate columns. Range metadata columns are written as regular
#' columns in the output.
#'
#' \strong{\code{GenomicRangesList} objects:} Separates list element metadata
#' from ranges data, writing them to different paths with schema properties
#' that describe how to split the ranges into list elements.
#'
#' \strong{\code{SelfHits} objects:} Writes graph edge lists as a data frame with
#' \code{from}, \code{to}, and optional metadata columns. The node count
#' (\code{nnode}) is stored in the schema properties to enable reconstruction as
#' a \linkS4class{DuckDBSelfHits} object.
#'
#' \strong{\code{SummarizedExperiment} objects:} Writes multi-assay experiments
#' with separate paths for feature data, sample data, and assay data.
#'
#' \strong{\code{RangedSummarizedExperiment} objects:} Extends
#' \code{SummarizedExperiment} functionality with genomic ranges for features.
#'
#' \strong{\code{SingleCellExperiment} objects:} Extends
#' \code{SummarizedExperiment} with additional single-cell specific data
#' including reduced dimensions, alternate experiments, and column/row pairings
#' (\code{colPairs}/\code{rowPairs}). Pairwise graphs are written to
#' \code{sample_graphs/} and \code{feature_graphs/} subdirectories.
#'
#' \strong{\code{MultiAssayExperiment} objects:} Writes multi-experiment studies
#' with separate paths for sample data, sample mapping, and experiment data.
#'
#' @section Frictionless Data Package Metadata:
#' This function implements the Frictionless Data Package specification for
#' metadata, creating \code{datapackage.json} files that describe the data
#' structures of the written objects. This function creates a separate
#' \code{datapackage.json} file for each object type, including:
#' \itemize{
#'   \item \code{name} - The name/class of the data package
#'   \item \code{resources} - A list of data resources, each containing:
#'   \itemize{
#'     \item \code{name} - Resource identifier
#'     \item \code{path} - Relative path to the data files
#'     \item \code{class} - Class of the resource (e.g., "data_frame", "array")
#'     \item \code{format} - File format ("parquet")
#'     \item \code{mediatype} - MIME type ("application/vnd.apache.parquet")
#'     \item \code{schema} - Field definitions with types and constraints
#'   }
#' }
#' This standardized metadata format ensures interoperability and provides
#' comprehensive documentation of the data structure for downstream analysis.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{readParquet}} for reading Bioconductor objects from parquet
#'   \item \code{\link{createDimTables}} for creating dimension lookup tables
#'   \item \code{\link[arrow]{write_dataset}} for the underlying Arrow functionality
#'   \item \code{\link[S4Arrays]{ArrayGrid}} for grid partitioning
#'   \item \code{\link[BiocParallel]{BiocParallelParam}} for parallel processing options
#' }
#'
#' @examples
#' # Write the Titanic dataset to a single Parquet file
#' tf1 <- tempfile()
#' writeParquet(Titanic, tf1)
#' list.files(tf1, full.names = TRUE, recursive = TRUE)
#'
#' # Write the state.x77 matrix to multiple Parquet files using grid partitioning
#' tf3 <- tempfile()
#' state_grid <- RegularArrayGrid(dim(state.x77), c(10, 4))
#' writeParquet(state.x77, file.path(tf3, "state"), grid = state_grid)
#' dimtbls <- createDimTables(state.x77, grid = state_grid)
#' list.files(tf3, full.names = TRUE, recursive = TRUE)
#'
#' @include fieldtypes.R
#'
#' @keywords IO
#'
#' @export
#' @import methods BiocGenerics
#' @rdname writeParquet
setGeneric("writeParquet", signature = "x",
function(x, path, ...)
{
  standardGeneric("writeParquet")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

.prefixSeq <- function(prefix, n) {
    sprintf(paste0(prefix, "%0", floor(log10(n)) + 1L, "d"), seq_len(n))
}

#' @importFrom S4Vectors metadata
.vectorMetadata <-
function(x)
{
    md <- metadata(x)
    md <- md[vapply(md, is.vector, logical(1L))]
    if (length(md) == 0L)
        md <- NULL
    md
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Array-like objects
###

#' @export
#' @importFrom DelayedArray blockApply currentViewport defaultAutoGrid
#' @importFrom DelayedArray effectiveGrid getAutoBPPARAM
#' @importFrom DuckDBArray writeCoordArray
#' @importFrom S4Arrays mapToGrid
#' @importFrom S4Vectors head tail
#' @importFrom SparseArray COO_SparseArray nzvals
#' @rdname writeParquet
setMethod("writeParquet", "ANY",
function(x,
         path,
         indexcols = names(dimnames(x)) %||% sprintf("index%d", seq_along(dim(x))),
         indexrefs = NULL,
         datacol = "value",
         grid = defaultAutoGrid(COO_SparseArray(dim(x))),
         grid_suffix = "_group",
         BPPARAM = getAutoBPPARAM(),
         name = basename(path),
         class = "array",
         ...)
{
    if (is.null(dim(x))) {
        stop("the default method of writeParquet requires 'x' to be array-like")
    }

    if (!(is.null(indexrefs) || length(indexrefs) == length(indexcols))) {
        stop("'indexrefs' must be NULL or a list of length(indexcols)")
    }

    if (inherits(x, "table")) {
        x <- unclass(x)
    }

    # Make column names unique
    unique_names <- make.unique(c(indexcols, datacol), sep = "_")
    indexcols <- head(unique_names, -1L)
    datacol <- tail(unique_names, 1L)

    # Manage dimnames
    dimnames_x <- dimnames(x) %||% lapply(dim(x), function(d) NULL)
    dimnames(x) <- lapply(dim(x), function(d) as.character(seq_len(d)))

    # Write array in coordinate format
    writeCoordArray(x, path = path, indexcols = indexcols, datacol = datacol,
                    grid = grid, grid_suffix = grid_suffix, BPPARAM = BPPARAM,
                    ...)

    # Generate field metadata
    fields <- c(lapply(seq_along(indexcols),
                       function(i) {
                           if (is.null(dimnames_x[[i]])) {
                               list(name = indexcols[i])
                           } else {
                               list(name = indexcols[i], type = "integer",
                                    categories =
                                        lapply(seq_along(dimnames_x[[i]]),
                                               function(j) {
                                                   list(value = j,
                                                        label = dimnames_x[[i]][[j]])
                                               }),
                                    categoriesOrdered = TRUE)
                           }
                       }),
                list(list(name = datacol)))

    # Generate foreign key metadata
    if (is.null(indexrefs)) {
        indexrefs <- rep(list(fields = ""), length(indexcols))
    }
    indexrefs <- lapply(seq_along(indexcols), function(i) {
        list(fields = indexcols[i], reference = indexrefs[[i]])
    })

    schema <- list(fields = fields, foreignKeys = indexrefs)
    resources <- list(list(name = name,
                           path = basename(path),
                           class = class,
                           format = "parquet",
                           mediatype = "application/vnd.apache.parquet",
                           schema = schema))

    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### data.frame & sf helper functions
###

.geometryCol <- function(x) {
    if (inherits(x, "sf")) {
        return(attr(x, "sf_column") %||% "geometry")
    }
    nms <- names(x)
    for (nm in c("geometry", "geom", "wkb")) {
        if (nm %in% nms) {
            col <- x[[nm]]
            if (inherits(col, "sfc"))
                return(nm)
            if (is.list(col) && length(col) > 0L &&
                all(vapply(col, function(r) is.null(r) || is.raw(r), NA)))
                return(nm)
        }
    }
    NULL
}

#' @importFrom S4Vectors I
.writeDataFrameParquet <-
function(x, path, indexcol, keycol, dimtbl, name, class, ...)
{
    if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

   if (is.null(indexcol)) {
        index <- NULL
    } else {
        index <- setNames(list(seq_len(nrow(x))), indexcol)
    }

    rnms <- attr(x, "row.names")
    if (is.null(keycol) || is.integer(rnms)) {
        key <- NULL
    } else {
        key <- setNames(list(rnms), keycol)
    }

    is_sf <- inherits(x, "sf")
    if (is_sf) {
        if (!requireNamespace("DuckDBSpatial", quietly = TRUE)) {
            stop("DuckDBSpatial package required for GeoParquet support; ",
                 "install with BiocManager::install('DuckDBSpatial')")
        }
        # Combine the columns into a single data.frame
        x <- do.call(cbind.data.frame, c(index, key, dimtbl, x))
        colnames(x) <- make.unique(colnames(x), sep = "_")
        x <- sf::st_sf(x)

        # defaults to compression = "zstd" and compression_level = 3L
        DuckDBSpatial::writeGeoParquet(x, file.path(path, "part-0.parquet"),
                         geom = attr(x, "sf_column"), ...)
    } else {
        # Protect the list columns from being coerced to atomic vectors
        for (j in seq_along(x)) {
            if (is.list(x[[j]])) {
                x[[j]] <- I(x[[j]])
            }
        }

        # Combine the columns into a single data.frame
        x <- do.call(cbind.data.frame, c(index, key, dimtbl, x))
        colnames(x) <- make.unique(colnames(x), sep = "_")

        # Optimize integer column storage
        for (j in seq_along(x)) {
            if (is.integer(x[[j]]) && length(x[[j]]) > 0L) {
                x[[j]] <- Array$create(x[[j]], type = .arrowType(x[[j]]))
            }
        }

        write_dataset(x, path, format = "parquet", compression = "zstd",
                      compression_level = 3L, ...)
    }

    schema <- list(fields = lapply(colnames(x), function(j) {
                       .buildFieldSpec(name = j, x = x[[j]])
                   }))

    # sortOrder represents physical organization (row index)
    if (!is.null(indexcol)) {
        schema[["sortOrder"]] <- list(
            list(field = indexcol, direction = "ascending")
        )
    }

    # primaryKey represents logical identity (row names if present, else row index)
    if (!is.null(key)) {
        schema[["primaryKey"]] <- keycol
    } else if (!is.null(indexcol)) {
        schema[["primaryKey"]] <- indexcol
    }

    # Add partitioning columns (if dimtbl provided)
    if (!is.null(dimtbl) && ncol(dimtbl) > 0L) {
        schema[["partitioning"]] <- colnames(dimtbl)
    }

    list(list(name = name,
              path = basename(path),
              class = class,
              format = "parquet",
              mediatype = "application/vnd.apache.parquet",
              schema = schema))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### data.frame objects
###

#' @export
#' @importFrom arrow Array write_dataset
#' @importFrom stats setNames
#' @rdname writeParquet
setMethod("writeParquet", "data.frame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "data.frame",
         ...)
{
    resources <- .writeDataFrameParquet(x, path = path, indexcol = indexcol,
                                        keycol = keycol, dimtbl = dimtbl,
                                        name = name, class = class, ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DataFrame objects
###

#' @export
#' @importClassesFrom S4Vectors DataFrame
#' @rdname writeParquet
setMethod("writeParquet", "DataFrame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "data_frame",
         ...)
{
    df <- as.data.frame(x, optional = TRUE)

    geom <- .geometryCol(df)
    is_sf <- !is.null(geom)
    if (is_sf) {
        df <- sf::st_as_sf(df, sf_column_name = geom)
    }

    resources <- .writeDataFrameParquet(df, path = path, indexcol = indexcol,
                                        keycol = keycol, dimtbl = dimtbl,
                                        name = name, class = class, ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TransposedDataFrame objects
###

#' @export
#' @importClassesFrom S4Vectors TransposedDataFrame
#' @rdname writeParquet
setMethod("writeParquet", "TransposedDataFrame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "transposed_data_frame",
         ...)
{
    resources <- callGeneric(t(x), path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, class = class, ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRanges objects
###

#' @export
#' @importClassesFrom GenomicRanges GenomicRanges
#' @rdname writeParquet
setMethod("writeParquet", "GenomicRanges",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "genomic_ranges",
         ...)
{
    # Write Data Table
    df <- as.data.frame(x, optional = TRUE)
    df[["width"]] <- NULL
    resources <- callGeneric(df, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, class = class, ...)

    # Add domain-specific metadata to schema
    schema <- resources[[length(resources)]][["schema"]]

    # Add genomic metadata and constraints
    schema <- .addGenomicMetadata(schema)

    resources[[length(resources)]][["schema"]] <- schema

    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRangesList objects
###

#' @export
#' @importClassesFrom GenomicRanges GenomicRangesList
#' @importClassesFrom IRanges CharacterList
#' @importFrom GenomicRanges seqnames start end strand
#' @importFrom S4Vectors DataFrame mcols
#' @rdname writeParquet
setMethod("writeParquet", "GenomicRangesList",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "genomic_ranges_list",
         ...)
{
    # Convert GenomicRangesList to DataFrame with AtomicList columns
    df <- DataFrame(seqnames = as(seqnames(x), "CharacterList"),
                    start = start(x), end = end(x),
                    strand = as(strand(x), "CharacterList"))
    rownames(df) <- names(x)

    if (is(x, "CompressedGenomicRangesList")) {
        # Add element-level metadata, if present
        elmcols <- mcols(unlist(x, use.names = FALSE))
        if (NCOL(elmcols) > 0L) {
            for (j in colnames(elmcols)) {
                df[[j]] <- relist(elmcols[[j]], x)
            }
        }
    }

    # Add list element-level metadata, if present
    if (NCOL(mcols(x)) > 0L) {
        df <- cbind(df, mcols(x))
        colnames(df) <- make.unique(colnames(df), sep = "_")
    }

    df <- as.data.frame(df, optional = TRUE)

    # Write Data Table
    resources <- callGeneric(df, path = path, indexcol = indexcol, keycol = keycol,
                             dimtbl = dimtbl, name = name, class = class, ...)

    # Add domain-specific metadata to schema for ranges
    schema <- resources[[length(resources)]][["schema"]]
    rnms <- sapply(schema[["fields"]], `[[`, "name")

    # Add genomic metadata and constraints
    schema <- .addGenomicMetadata(schema)

    resources[[length(resources)]][["schema"]] <- schema

    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SelfHits objects
###

#' @export
#' @importClassesFrom S4Vectors SelfHits
#' @importFrom S4Vectors nnode
#' @rdname writeParquet
setMethod("writeParquet", "SelfHits",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         name = basename(path),
         class = "graph_edges",
         ...)
{
    df <- as.data.frame(x, optional = TRUE)
    resources <- callGeneric(df, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, class = class, ...)

    # Add domain-specific metadata to schema
    schema <- resources[[length(resources)]][["schema"]]

    # Add graph metadata and constraints
    schema <- .addGraphMetadata(schema, nnode(x))

    resources[[length(resources)]][["schema"]] <- schema

    invisible(resources)
})

.writeParquetGraphs <-
function(x,
         path,
         package = list(class = "graphs", resources = list()),
         ...)
{
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]

        if (is(x_i, "DualSubset")) {
            x_i <- x_i@hits
        }

        path_i <- paste0("graph=", nms_i)
        resources <- writeParquet(x_i, path = file.path(path, path_i),
                                  name = nms_i, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assays objects
###

#' @export
#' @importClassesFrom S4Vectors DataFrame TransposedDataFrame
#' @importClassesFrom SummarizedExperiment Assays
#' @importFrom S4Vectors getListElement
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "Assays",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         indexrefs = NULL,
         datacol = "value",
         grid = defaultAutoGrid(COO_SparseArray(dim(x[[1L]]))),
         grid_suffix = "group__",
         package = list(class = "assays", resources = list()),
         ...)
{
    # Write transpose of the assays
    indexcols <- rev(indexcols)
    indexrefs <- rev(indexrefs)
    grid <- t(grid)
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- t(getListElement(x, i))
        path_i <- paste0("assay=", nms_i)
        if (is(x_i, "DataFrame")) {
            resources <-
                callGeneric(x_i, path = file.path(path, path_i),
                            indexcol = indexcols[1L], name = nms_i, ...)
        } else if (is(x_i, "TransposedDataFrame")) {
            resources <-
                callGeneric(x_i, path = file.path(path, path_i),
                            indexcol = indexcols[2L], name = nms_i, ...)
        } else {
            # Array-like object
            dimnames(x_i) <- NULL
            resources <-
                callGeneric(x_i, path = file.path(path, path_i),
                            indexcols = indexcols, indexrefs = indexrefs,
                            datacol = datacol, grid = grid,
                            grid_suffix = grid_suffix, name = nms_i, ...)
        }
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SummarizedExperiment objects
###

#' @export
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom DelayedArray defaultAutoGrid
#' @importFrom DuckDBArray createDimTables
#' @importFrom SparseArray COO_SparseArray
#' @importFrom SummarizedExperiment assays colData rowData rowRanges
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "SummarizedExperiment",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = ifelse(is(x, "RangedSummarizedExperiment"),
                                       "ranged_summarized_experiment",
                                       "summarized_experiment"),
                        resources = list()),
         ...)
{
    # Make dimnames unique
    rownames(x) <- make.unique(rownames(x), sep = "_")
    colnames(x) <- make.unique(colnames(x), sep = "_")

    # Dimension Tables
    grid <- defaultAutoGrid(COO_SparseArray(dim(assays(x)[[1L]])))
    grid_suffix <- "group__"
    dimtbls <- createDimTables(assays(x)[[1L]],
                               indexcols = indexcols,
                               grid = grid,
                               grid_suffix = grid_suffix)

    # Feature Data
    if (is(x, "RangedSummarizedExperiment")) {
        features <- rowRanges(x)
    } else {
        features <- rowData(x)
    }
    resources <-
        callGeneric(features, path = file.path(path, "features"),
                    indexcol = indexcols[1L], dimtbl = dimtbls[[1L]], ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Sample Data
    resources <-
        callGeneric(colData(x), path = file.path(path, "samples"),
                    indexcol = indexcols[2L], dimtbl = dimtbls[[2L]], ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Assay Data
    resources <- list(list(name = "assays",
                           path = "assays",
                           class = "data_package"))
    indexrefs <- list(list(fields = indexcols[1L],
                           reference = list(fields = indexcols[1L],
                                            resource = "../features")),
                      list(fields = indexcols[2L],
                           reference = list(fields = indexcols[2L],
                                            resource = "../samples")))
    callGeneric(x@assays, path = file.path(path, resources[[1L]][["path"]]),
                indexcols = indexcols, indexrefs = indexrefs, grid = grid,
                grid_suffix = grid_suffix, ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    package[["annotations"]] <- .vectorMetadata(x)

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Modalities objects
###

#' @importFrom jsonlite write_json
.writeParquetModalities <-
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = "modalities", resources = list()),
         ...)
{
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]
        path_i <- paste0("modality=", nms_i)
        resources <- list(name = nms_i, path = path_i, class = "data_package")
        resources <- list(resources)
        writeParquet(x_i, path = file.path(path, path_i), indexcols = indexcols,
                     ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SingleCellExperiment objects
###

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors I
#' @importFrom SingleCellExperiment altExps mainExpName
#' @importFrom SingleCellExperiment reducedDims reducedDimNames
#' @importFrom SingleCellExperiment colPairs rowPairs
#' @rdname writeParquet
setMethod("writeParquet", "SingleCellExperiment",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = "single_cell_experiment",
                        resources = list()),
         ...)
{
    # Make dimnames unique
    rownames(x) <- make.unique(rownames(x), sep = "_")
    colnames(x) <- make.unique(colnames(x), sep = "_")

    # Row Loadings
    loadings <- rowLoadings(x)
    if (length(loadings)) {
        loadings <- lapply(seq_along(loadings), function(i) {
            mat <- as.matrix(loadings[[i]])
            dimnames(mat) <- NULL
            I(asplit(mat, 1L))
        })
        loadings <- do.call(data.frame, loadings)
        rownames(loadings) <- rownames(x)
        colnames(loadings) <- make.unique(rowLoadingNames(x), sep = "_")
        resources <- callGeneric(loadings,
                                 path = file.path(path, "feature_embeddings"),
                                 indexcol = indexcols[1L], class = "data_frame",
                                 ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Reduced Dimensions
    rdims <- reducedDims(x)
    if (length(rdims)) {
        rdims <- lapply(seq_along(rdims), function(j) {
            mat <- as.matrix(rdims[[j]])
            dimnames(mat) <- NULL
            I(asplit(mat, 1L))
        })
        rdims <- do.call(data.frame, rdims)
        rownames(rdims) <- colnames(x)
        colnames(rdims) <- make.unique(reducedDimNames(x), sep = "_")
        resources <- callGeneric(rdims,
                                 path = file.path(path, "sample_embeddings"),
                                 indexcol = indexcols[2L], class = "data_frame",
                                 ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Alternative Experiments
    exps <- altExps(x)
    if (length(exps)) {
        resources <- list(list(name = "modalities",
                               path = "modalities",
                               class = "data_package"))
        .writeParquetModalities(exps,
                                path = file.path(path, resources[[1L]][["path"]]),
                                indexcols = indexcols, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Row Pairs
    rpairs <- rowPairs(x, asSparse = FALSE)
    if (length(rpairs)) {
        resources <- list(list(name = "feature_graphs",
                               path = "feature_graphs",
                               class = "data_package"))
        .writeParquetGraphs(rpairs,
                            path = file.path(path, resources[[1L]][["path"]]),
                            ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Column Pairs
    cpairs <- colPairs(x, asSparse = FALSE)
    if (length(cpairs)) {
        resources <- list(list(name = "sample_graphs",
                               path = "sample_graphs",
                               class = "data_package"))
        .writeParquetGraphs(cpairs,
                            path = file.path(path, resources[[1L]][["path"]]),
                            ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Main Experiment Name
    package[["main_exp_name"]] <- mainExpName(x)

    callNextMethod(x, path = path, indexcols = indexcols, package = package, ...)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ExperimentList objects
###

#' @export
#' @importClassesFrom MultiAssayExperiment ExperimentList
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "ExperimentList",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = "experiment_list", resources = list()),
         ...)
{
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]
        path_i <- paste0("experiment=", nms_i)
        if (is(x_i, "SummarizedExperiment")) {
            resources <- list(name = nms_i, path = path_i, class = "data_package")
            resources <- list(resources)
            callGeneric(x_i, path = file.path(path, path_i), indexcols = indexcols, ...)
        } else {
            # Array-like object
            resources <-
                writeParquet(x_i, path = file.path(path, path_i),
                             indexcols = indexcols, datacol = "value",
                             grid_suffix = "group__", name = nms_i, ...)
        }
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MultiAssayExperiment objects
###

#' @export
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment experiments sampleMap
#' @importFrom SummarizedExperiment colData
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "MultiAssayExperiment",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = "multi_assay_experiment",
                        resources = list()),
         ...)
{
    # Subject Data
    resources <- callGeneric(colData(x), path = file.path(path, "subjects"),
                             ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Sample Map
    resources <- callGeneric(sampleMap(x), path = file.path(path, "sample_map"),
                             ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Experiment Data
    resources <- list(list(name = "experiments",
                           path = "experiments",
                           class = "experiment_list"))
    callGeneric(experiments(x), path = file.path(path, resources[[1L]][["path"]]),
                indexcols = indexcols, ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    package[["annotations"]] <- .vectorMetadata(x)

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PointsLayerList objects
###

#' @export
#' @importClassesFrom MultiAssaySpatialExperiment PointsLayerList
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "PointsLayerList",
function(x, path, name = basename(path),
         class = "points_layer_list",
         package = list(class = "points_layer_list", resources = list()),
         ...)
{
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]
        path_i <- paste0("point=", nms_i)
        resources <- callGeneric(x_i, path = file.path(path, path_i), name = nms_i, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ShapesLayerList objects
###

#' @export
#' @importClassesFrom MultiAssaySpatialExperiment ShapesLayerList
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "ShapesLayerList",
function(x,
         path,
         name = basename(path),
         class = "shapes_layer_list",
         package = list(class = "shapes_layer_list", resources = list()),
         ...)
{
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]
        path_i <- paste0("shape=", nms_i)
        resources <- callGeneric(x_i, path = file.path(path, path_i), name = nms_i, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MultiAssaySpatialExperiment objects
###

#' @export
#' @importClassesFrom MultiAssaySpatialExperiment MultiAssaySpatialExperiment
#' @importFrom jsonlite write_json
#' @rdname writeParquet
setMethod("writeParquet", "MultiAssaySpatialExperiment",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(class = "multi_assay_spatial_experiment",
                        resources = list()),
         ...)
{
    # Points
    pts <- MultiAssaySpatialExperiment::spatialPoints(x)
    if (length(pts) > 0L) {
        resources <- list(list(name = "points",
                               path = "points",
                               class = "points_layer_list"))
        callGeneric(pts, path = file.path(path, resources[[1L]][["path"]]), ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Shapes
    shps <- MultiAssaySpatialExperiment::spatialShapes(x)
    if (length(shps) > 0L) {
        resources <- list(list(name = "shapes",
                               path = "shapes",
                               class = "shapes_layer_list"))
        callGeneric(shps, path = file.path(path, resources[[1L]][["path"]]), ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Images path references
    imgs <- MultiAssaySpatialExperiment::spatialImages(x)
    if (length(imgs) > 0L) {
        imgs_dir <- file.path(path, "images")
        dir.create(imgs_dir, showWarnings = FALSE, recursive = TRUE)
        for (i in seq_along(imgs)) {
            nm <- names(imgs)[i]
            img_meta <- list(name = nm, type = "image", path = NA_character_)
            write_json(img_meta, file.path(imgs_dir, paste0(nm, ".json")),
                       auto_unbox = TRUE, pretty = TRUE)
        }
        package[["resources"]] <- c(package[["resources"]],
                                    list(list(name = "images",
                                              path = "images",
                                              class = "raster_layer_list")))
    }

    # Labels path references
    lbls <- MultiAssaySpatialExperiment::spatialLabels(x)
    if (length(lbls) > 0L) {
        lbls_dir <- file.path(path, "labels")
        dir.create(lbls_dir, showWarnings = FALSE, recursive = TRUE)
        for (i in seq_along(lbls)) {
            nm <- names(lbls)[i]
            lbl_meta <- list(name = nm, type = "label", path = NA_character_)
            write_json(lbl_meta, file.path(lbls_dir, paste0(nm, ".json")),
                       auto_unbox = TRUE, pretty = TRUE)
        }
        package[["resources"]] <- c(package[["resources"]],
                                    list(list(name = "labels",
                                              path = "labels",
                                              class = "raster_layer_list")))
    }

    # Image Data
    img_data <- MultiAssaySpatialExperiment::imgData(x)
    if (!is.null(img_data) && nrow(img_data) > 0L) {
        # Create images directory for materialized image files
        imgs_dir <- file.path(path, "imgdata_images")
        dir.create(imgs_dir, showWarnings = FALSE, recursive = TRUE)

        if ("data" %in% colnames(img_data)) {
            img_data_copy <- img_data
            img_objs <- img_data[["data"]]

            # Save each image and create relative path references
            img_paths <- vapply(seq_along(img_objs), function(i) {
                obj <- img_objs[[i]]
                if (is.null(obj)) {
                    return(NA_character_)
                }

                # Filename: imgdata_images/img_<i>.png
                img_filename <- sprintf("img_%d.png", i - 1L)
                img_filepath <- file.path(imgs_dir, img_filename)

                # Try to copy existing file if it's a StoredSpatialImage
                if (is(obj, "StoredSpatialImage")) {
                    src <- tryCatch(
                        SpatialExperiment::imgSource(obj),
                        error = function(e) NULL
                    )
                    if (!is.null(src) && file.exists(src)) {
                        # Check if it's already PNG, otherwise convert
                        if (grepl("\\.png$", src, ignore.case = TRUE)) {
                            file.copy(src, img_filepath, overwrite = TRUE)
                            return(file.path("imgdata_images", img_filename))
                        }
                    }
                }

                # Fallback: extract raster and save as PNG
                tryCatch({
                    if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
                        ras <- SpatialExperiment::imgRaster(obj)
                        # Convert raster to PNG format
                        Y <- grDevices::col2rgb(as.matrix(ras))
                        Y <- t(Y)
                        Y <- Y / 255
                        dim(Y) <- c(dim(ras), ncol(Y))

                        if (requireNamespace("png", quietly = TRUE)) {
                            png::writePNG(Y, target = img_filepath)
                            file.path("imgdata_images", img_filename)
                        } else {
                            warning("png package not available; skipping image ", i)
                            NA_character_
                        }
                    } else {
                        NA_character_
                    }
                }, error = function(e) {
                    warning("Failed to save image ", i, ": ", e$message)
                    NA_character_
                })
            }, character(1L))

            # Store relative paths to materialized images
            img_data_copy[["image_file"]] <- img_paths

            # Remove original data column (not serializable)
            img_data_copy[["data"]] <- NULL

            resources_img <- callGeneric(img_data_copy, path = file.path(path, "img_data"),
                                         ...)
        } else {
            resources_img <- callGeneric(img_data, path = file.path(path, "img_data"),
                                         ...)
        }
        package[["resources"]] <- c(package[["resources"]], resources_img)
    }

    # Spatial Map
    spatial_map <- MultiAssaySpatialExperiment::spatialMap(x)
    if (!is.null(spatial_map) && nrow(spatial_map) > 0L) {
        resources_sm <- callGeneric(spatial_map, path = file.path(path, "spatial_map"),
                                    ...)
        package[["resources"]] <- c(package[["resources"]], resources_sm)
    }

    callNextMethod(x, path = path, indexcols = indexcols, package = package, ...)

    invisible(NULL)
})
