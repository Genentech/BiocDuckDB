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
#'   \item \code{data.frame} / \code{DataFrame} objects - written with field type metadata
#'   \item \code{DuckDBTable} / \code{DuckDBDataFrame} objects - lazy SQL
#'     \code{COPY TO} export via \code{\link[DuckDBDataFrame]{writeDuckDBTableParquet}}
#'     in \code{\link[DuckDBDataFrame]{parquet-io}} (same \code{append}, \code{offset},
#'     \code{part} contract as in-memory tables)
#'   \item \code{DuckDBTransposedDataFrame}, \code{DuckDBSelfHits},
#'     \code{DuckDBGRanges}, \code{DuckDBGRangesList} - delegate to the underlying
#'     lazy table without materializing
#'   \item \code{TransposedDataFrame} objects - transposed before writing
#'   \item \code{list} / \code{List} objects - each element dispatched individually;
#'     requires \code{dimension}; \code{layout} defaults to \code{NULL} (each
#'     element's method supplies its own default layout)
#'   \item \code{GenomicRanges} objects - genomic coordinates with metadata
#'   \item \code{GenomicRangesList} objects - lists of genomic ranges stored as
#'     LIST-typed columns
#'   \item \code{SelfHits} objects - graph edge lists with node count metadata
#'   \item \code{Assays} objects - collection of feature-by-sample matrices
#'   \item \code{SummarizedExperiment} objects - feature/sample metadata with assays
#'   \item \code{RangedSummarizedExperiment} objects - as above with genomic ranges
#'     for features
#'   \item \code{SingleCellExperiment} objects - extends \code{SummarizedExperiment}
#'     with embeddings, alt experiments, dimension tables, and pairwise graphs
#'   \item \code{ExperimentList} objects - named collection of experiments;
#'     primitive method used internally by \code{MultiAssayExperiment} and
#'     \code{SingleCellExperiment} (altExps); returns resources without writing
#'     \code{datapackage.json}
#'   \item \code{MultiAssayExperiment} objects - multi-experiment study with subject
#'     metadata and sample map
#'   \item \code{PointsLayerList} / \code{ShapesLayerList} objects - spatial point
#'     and polygon layers
#'   \item \code{MultiAssaySpatialExperiment} objects - spatially-resolved
#'     multi-assay experiment
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
#' @param refs For flat (\code{data.frame} / \code{DataFrame}) tables, an
#' optional list of foreign key references to emit in the schema, the flat
#' analog of \code{indexrefs}. Each entry is
#' \code{list(fields = <local>, reference = list(fields = <target>, resource = <res>))}.
#' Used, for example, to declare the \code{MultiAssayExperiment} sample_map
#' \code{primary} -> \code{subjects} correspondence. Defaults to \code{NULL}.
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
#' @param dimension A character string describing which biological axis the
#' resource is indexed by. One of \code{"feature"}, \code{"sample"},
#' \code{"crossed"} (feature \eqn{\times} sample Cartesian product), or
#' \code{"unbound"} (not tied to either axis).
#' @param layout A character string identifying the physical storage layout of
#' the resource. Recorded in \code{datapackage.json} and used by
#' \code{readParquet} to dispatch the correct reader. Examples:
#' \code{"data_frame"}, \code{"coord_array"}, \code{"embedding_table"},
#' \code{"genomic_ranges"}, \code{"graph_edges"}, \code{"nested_data_frame"},
#' \code{"nested_experiment"}.
#' @param append Logical. For \code{data.frame}, \code{DataFrame}, and
#' \code{DuckDBTable} methods: when \code{TRUE}, append a new flat
#' \code{part-*.parquet} file under \code{path} (requires \code{part}). For
#' array-like objects (\code{ANY} method): forwarded to
#' \code{\link[DuckDBArray]{writeCoordArray}} for hive-partition append
#' (requires \code{length(grid) > 1L}; see \code{along}, \code{group_offset}
#' there). Defaults to \code{FALSE}.
#' @param offset Non-negative integer. For \code{data.frame} and
#' \code{DuckDBTable} methods: added to the index column when \code{indexcol}
#' is set (\code{offset + seq_len(nrow(x))} in memory;
#' \code{offset + row_number()} in SQL for lazy tables); ignored when
#' \code{indexcol} is \code{NULL}. For array-like objects: forwarded to
#' \code{\link[DuckDBArray]{writeCoordArray}} as the coordinate shift along
#' \code{along}. Defaults to \code{0L}.
#' @param part Optional non-negative integer; \code{data.frame} and
#' \code{DuckDBTable} methods. When set, write a single \code{part-<n>.parquet}
#' file (via \code{arrow::write_parquet} for in-memory tables, or DuckDB
#' \code{COPY TO} for lazy tables). Required when \code{append = TRUE} on table
#' writes.
#' @param part_digits Zero-padding width for \code{part} in the filename (e.g.
#' \code{2L} yields \code{part-00.parquet}). Table methods (in-memory and lazy).
#' @param package For experiment-level methods (\code{SummarizedExperiment},
#' etc.), a list used to accumulate \code{datapackage.json} contents before
#' writing. Primitive \code{data.frame} methods ignore this argument; capture
#' the return value on the first flat part and append later parts with
#' \code{append = TRUE}. Contains:
#' \itemize{
#'   \item \code{model} - Package-level schema identifier (e.g.
#'     \code{"single_cell_experiment"}). Determines which \code{readParquet}
#'     reader reconstructs the object. Distinct from the resource-level
#'     \code{layout} field — \code{model} describes the whole package;
#'     \code{layout} describes one resource within it.
#'   \item \code{resources} - Accumulating list of resource entries, each
#'     contributed by a primitive \code{writeParquet} method call.
#' }
#' Callers should not normally set this; the default is supplied by each
#' experiment method. Experiment methods accumulate \code{package[["resources"]]}
#' within a single call and write \code{datapackage.json} at the end. Primitive
#' \code{data.frame} methods do not update a caller's \code{package} across
#' separate invocations; capture the list returned from the first
#' \code{writeParquet} call and merge manually (see \strong{Flat table append}
#' below).
#' @param ... Additional arguments. For \code{data.frame} methods with
#' \code{part} set: passed to \code{\link[arrow]{write_parquet}}. Otherwise
#' for tables: passed to \code{\link[arrow]{write_dataset}}. For array-like
#' objects (\code{ANY} method): passed to
#' \code{\link[DuckDBArray]{writeCoordArray}}, including:
#' \describe{
#'   \item{\code{append}}{Logical; hive-partition append when \code{TRUE}
#'     (requires \code{length(grid) > 1L}). See \code{?writeCoordArray}.}
#'   \item{\code{along}}{Integer; append dimension (required when
#'     \code{append = TRUE}).}
#'   \item{\code{group_offset}}{Non-negative integer; partition-group offset
#'     along \code{along}; defaults to \code{0L}.}
#'   \item{\code{arrowtype}}{Optional \code{\link[arrow]{DataType}} for the
#'     value column; schema pinning on append is described in
#'     \code{?writeCoordArray}.}
#'   \item{\code{max_dim}}{Optional length-\code{ndim(x)} integer vector of
#'     index upper bounds.}
#'   \item{\code{existing_data_behavior}}{Passed through to
#'     \code{\link[arrow]{write_dataset}} when applicable.}
#' }
#' Note: \code{offset} for array-like objects is documented above; table
#' methods use \code{offset} as an explicit formal argument.
#'
#' @return
#' For primitive methods (\code{ANY}, \code{data.frame}, \code{DataFrame},
#' \code{TransposedDataFrame}, \code{GenomicRanges}, \code{GenomicRangesList},
#' \code{SelfHits}, \code{Assays}, \code{list}, \code{List},
#' \code{ExperimentList}, \code{PointsLayerList}, \code{ShapesLayerList}):
#' invisibly returns a list of Frictionless resource entries suitable for
#' inclusion in a \code{datapackage.json} \code{resources} array.
#'
#' For experiment methods (\code{SummarizedExperiment},
#' \code{SingleCellExperiment}, \code{MultiAssayExperiment},
#' \code{MultiAssaySpatialExperiment}): invisibly returns \code{NULL}; the
#' \code{datapackage.json} is written to \code{path} as a side effect.
#'
#' For flat multi-part table writes (\code{data.frame} with \code{append = TRUE}
#' and \code{part > 0}): invisibly returns \code{NULL}.
#'
#' @details
#' This function provides specialized handling for different object types:
#'
#' \strong{Array-like objects:} Converts multi-dimensional arrays into a
#' coordinate (long) format where each non-zero element is represented as a
#' row with columns for each dimension and the value. For sparse arrays, only
#' non-zero elements are written, making it efficient for sparse data. When a
#' grid is provided with multiple cells, the array is partitioned and each
#' partition is written to a separate subdirectory. Array writes delegate to
#' \code{\link[DuckDBArray]{writeCoordArray}}; arguments
#' \code{append}, \code{along}, \code{offset}, \code{group_offset},
#' \code{arrowtype}, and \code{max_dim} are forwarded via \code{...} and
#' documented on the \code{writeCoordArray} help page.
#'
#' \strong{Flat table append (\code{data.frame} / \code{DataFrame}):} Use
#' \code{part}, \code{part_digits}, \code{append}, and \code{offset} to stream
#' chunked sample or feature tables as \code{part-0.parquet}, \code{part-1.parquet},
#' \ldots without temporary directories. The global index column (when
#' \code{indexcol} is set) uses \code{offset + seq_len(nrow(x))} on each chunk.
#' The first call returns a Frictionless resource list (one element); later
#' append parts return \code{NULL}. Build \code{datapackage.json} from the first
#' return, for example:
#' \preformatted{
#' res <- writeParquet(chunk1, samples_dir, indexcol = "__sample__",
#'                     part = 0L, name = "samples", dimension = "sample")
#' pkg <- list(resources = res)
#' writeParquet(chunk2, samples_dir, indexcol = "__sample__",
#'              offset = nrow(chunk1), part = 1L, append = TRUE, ...)
#' }
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
#' \strong{\code{SelfHits} objects:} Writes graph edge lists as a data frame
#' with \code{from}, \code{to}, and optional metadata columns. The node count
#' (\code{nnode}) is stored in the schema properties to enable reconstruction as
#' a \linkS4class{DuckDBSelfHits} object.
#'
#' \strong{Automatic unnesting of nested DataFrames:} When writing
#' \code{SingleCellExperiment} objects, any DataFrame-class columns found in
#' \code{rowData()} or \code{colData()} are automatically extracted and moved to
#' \code{rowTables()} or \code{colTables()} respectively. This unnesting enables
#' independent Parquet serialization of multi-valued properties (e.g., multiple
#' diseases per patient, multiple isoforms per gene) while keeping the main
#' metadata tables flat and SQL-queryable. The original nested columns are
#' removed from \code{rowData()}/\code{colData()} after extraction. This
#' behavior is automatic and requires no user intervention.
#'
#' \strong{\code{SummarizedExperiment} objects:} Writes multi-assay experiments
#' with separate paths for feature data, sample data, and assay data.
#'
#' \strong{\code{RangedSummarizedExperiment} objects:} Extends
#' \code{SummarizedExperiment} functionality with genomic ranges for features.
#'
#' \strong{\code{SingleCellExperiment} objects:} Extends
#' \code{SummarizedExperiment} with single-cell specific data. All collections
#' are written as flat resource directories directly under \code{path}. The
#' directory name encodes both axis and type: \code{feature_embeddings/},
#' \code{sample_embeddings/}, \code{feature_table_<name>/},
#' \code{sample_table_<name>/}, \code{feature_graph_<name>/},
#' \code{sample_graph_<name>/}, \code{experiment_<name>/}. The path prefix is
#' derived automatically from \code{layout}. Alternative experiments are written
#' directly to root (flattened, same pattern as \code{MultiAssayExperiment}
#' experiments). See the automatic unnesting section for details about
#' \code{DataFrame} columns in \code{rowData()}/\code{colData()}.
#'
#' \strong{\code{MultiAssayExperiment} objects:} Writes multi-experiment studies
#' with separate paths for sample data, sample mapping, and experiment data.
#'
#' @section Frictionless Data Package Metadata:
#' This function writes a \code{datapackage.json} file (Frictionless Data
#' Package v2.0) alongside the Parquet files. The file contains two levels of
#' metadata:
#'
#' \strong{Package-level fields} (top of \code{datapackage.json}):
#' \itemize{
#'   \item \code{model} - Overall schema identifier. Determines which
#'     \code{readParquet} reader is used to reconstruct the container object
#'     (e.g. \code{"single_cell_experiment"}, \code{"multi_assay_experiment"}).
#'     When absent, \code{readParquet} returns a \code{SimpleList} of resources
#'     with no imposed schema.
#'   \item \code{annotations} - Non-relational elements from \code{metadata(x)}
#'     (scalars, vectors, 1-D arrays, JSON-safe nested lists). Tabular /
#'     relational items (\code{DataFrame}, \code{matrix}, 2-D arrays) are
#'     written as \code{unbound} Parquet sidecars and referenced here via
#'     \code{parquet_ref} stubs; mixed nested lists use a
#'     \code{nested_mapping} wrapper.
#' }
#'
#' \strong{Per-resource fields} (each entry in the \code{resources} array):
#' \itemize{
#'   \item \code{name} - Resource identifier (typically the object name)
#'   \item \code{path} - Relative directory path to the Parquet dataset
#'   \item \code{dimension} - Biological axis: \code{"feature"}, \code{"sample"},
#'     \code{"crossed"}, or \code{"unbound"}
#'   \item \code{layout} - Physical storage layout (BiocDuckDB extension);
#'     together with \code{dimension}, uniquely determines the R accessor and
#'     AnnData slot used during reconstruction
#'   \item \code{format} - File format (\code{"parquet"})
#'   \item \code{mediatype} - MIME type (\code{"application/vnd.apache.parquet"})
#'   \item \code{schema} - Field definitions including types, primary key,
#'     sort order, and foreign key references
#' }
#' See the BiocDuckDB storage patterns documentation for the full
#' \code{model} table and \code{dimension} + \code{layout} dispatch table.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{readParquet}} for reading Bioconductor objects from parquet
#'   \item \code{\link[DuckDBArray]{writeCoordArray}} for coord-array layout,
#'     hive partitioning, and array append (\code{append}, \code{along},
#'     \code{offset}, \code{group_offset})
#'   \item \code{\link[DuckDBDataFrame]{parquet-io}} for shared Parquet I/O:
#'     flat append validation (\code{setupFlatParquetWrite}, \code{checkAppendPart},
#'     \code{validateAppendOffset}), SQL \code{COPY TO} helpers
#'     (\code{buildParquetCopySQL}), and lazy table export
#'     (\code{writeDuckDBTableParquet})
#'   \item \code{\link{createDimTables}} for creating dimension lookup tables
#'   \item \code{\link[arrow]{write_dataset}} and \code{\link[arrow]{write_parquet}}
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

### Frictionless profile that the emitted datapackage.json declares conformance
### to via its top-level "$schema". Pre-release we point at the canonical,
### resolvable base Data Package v2 profile: BiocDuckDB's extensions are
### additive and base-conformant by design, so this is a true, dereferenceable
### claim. Switch to the published BiocDuckDB profile $id at Bioconductor release
### (a one-line change here).
.BIOCDUCKDB_PROFILE <- "https://datapackage.org/profiles/2.0/datapackage.json"

.prefixSeq <- function(prefix, n) {
    sprintf(paste0(prefix, "%0", floor(log10(n)) + 1L, "d"), seq_len(n))
}

### Metadata serialization policy (relational vs non-relational):
###   JSON  — scalar-like values and 1-D sequences that jsonlite can encode,
###           plus plain nested lists whose children are all JSON-safe.
###   Parquet — data.frame, DFrame, matrix, arrays with dim >= 2, and other
###           types supported by writeParquet methods.
###   Skip  — unsupported S4 objects (warning issued).
###
### R's type predicates are a poor match for this policy: Date, POSIXt, and
### factor are not is.atomic(); package_version is is.list() but is a
### list-shaped scalar (length-1, self-referential), not a nested mapping.
### We therefore classify by structure (tabular / plain list / else) and
### use jsonlite to decide whether the remainder is a JSON leaf.

#' @importClassesFrom S4Vectors DFrame List
.isMetadataTabular <- function(x) {
    is.data.frame(x) || is(x, "DFrame") || is.matrix(x) ||
        (is.array(x) && !is.null(dim(x)) && length(dim(x)) >= 2L)
}

.isListShapedScalar <- function(x) {
    is.list(x) && length(x) == 1L && identical(x[[1L]], x)
}

.isPlainMetadataList <- function(x) {
    (is.list(x) || is(x, "List")) && !is.data.frame(x) && !is(x, "DFrame") &&
        !isS4(x) && !.isListShapedScalar(x)
}

.metadataJsonValue <- function(x) {
    if (is.null(x))
        return(NULL)
    if (is.array(x) && length(dim(x)) <= 1L)
        return(as.vector(x))
    if (is.atomic(x) && is.null(dim(x)) && !is.object(x))
        return(x)
    if (.isListShapedScalar(x))
        return(as.character(x))
    if (is.object(x) && !isS4(x) && !.isPlainMetadataList(x)) {
        ch <- tryCatch(as.character(x), error = function(e) NULL)
        if (!is.null(ch))
            return(ch)
    }
    x
}

.isMetadataJsonLeaf <- function(x) {
    if (is.null(x) || .isMetadataTabular(x) || .isPlainMetadataList(x) || isS4(x))
        return(is.null(x))
    tryCatch({
        jsonlite::toJSON(.metadataJsonValue(x), auto_unbox = TRUE, null = "null")
        TRUE
    }, error = function(e) FALSE)
}

.metadataTabularValue <- function(x) {
    if (is.matrix(x)) {
        df <- as.data.frame(x, stringsAsFactors = FALSE)
        if (!is.null(colnames(x)))
            colnames(df) <- colnames(x)
        df
    } else if (is.array(x) && length(dim(x)) >= 2L) {
        .metadataTabularValue(as.matrix(x))
    } else {
        x
    }
}

#' @importFrom S4Vectors metadata
.serializeMetadataValue <- function(x, path, name, resources, ...) {
    if (.isMetadataJsonLeaf(x)) {
        return(list(value = .metadataJsonValue(x), resources = resources,
                    has_parquet = FALSE))
    }
    if (.isMetadataTabular(x)) {
        path_k <- file.path(path, sprintf("unbound_%s", name))
        res <- try(writeParquet(.metadataTabularValue(x), path = path_k,
                                name = name, dimension = "unbound",
                                layout = "data_frame", ...),
                   silent = TRUE)
        if (inherits(res, "try-error")) {
            warning("Skipping metadata item '", name, "': ", res,
                    call. = FALSE)
            return(list(value = NULL, resources = resources,
                        has_parquet = FALSE, skip = TRUE))
        }
        stub <- list("__type__" = "parquet_ref", resource = name)
        return(list(value = stub, resources = c(resources, res),
                    has_parquet = TRUE))
    }
    if (.isPlainMetadataList(x)) {
        children <- list()
        child_has_parquet <- FALSE
        for (k in seq_along(x)) {
            sub_nm <- names(x)[k] %||% sprintf("item-%d", k)
            qualified <- sprintf("%s__%s", name, sub_nm)
            elt <- x[[k]]
            child <- if (identical(elt, x) || .isListShapedScalar(elt)) {
                list(value = .metadataJsonValue(x), resources = resources,
                     has_parquet = FALSE)
            } else {
                .serializeMetadataValue(elt, path = path, name = qualified,
                                        resources = resources, ...)
            }
            if (!isTRUE(child$skip)) {
                children[[sub_nm]] <- child$value
            }
            resources <- child$resources
            child_has_parquet <- child_has_parquet || isTRUE(child$has_parquet)
        }
        if (length(children) == 0L) {
            return(list(value = NULL, resources = resources,
                        has_parquet = FALSE, skip = TRUE))
        }
        value <- if (child_has_parquet) {
            c(list("__type__" = "nested_mapping"), children)
        } else {
            children
        }
        return(list(value = value, resources = resources,
                    has_parquet = child_has_parquet))
    }
    path_k <- file.path(path, sprintf("unbound_%s", name))
    res <- try(writeParquet(x, path = path_k, name = name,
                            dimension = "unbound", ...),
               silent = TRUE)
    if (inherits(res, "try-error")) {
        warning("Skipping unsupported metadata item '", name, "': ", res,
                call. = FALSE)
        return(list(value = NULL, resources = resources,
                    has_parquet = FALSE, skip = TRUE))
    }
    stub <- list("__type__" = "parquet_ref", resource = name)
    list(value = stub, resources = c(resources, res), has_parquet = TRUE)
}

#' @importFrom S4Vectors metadata
.serializeMetadata <- function(x, path, resources = list(), ...) {
    md <- metadata(x)
    if (length(md) == 0L)
        return(list(annotations = NULL, resources = resources))
    annotations <- list()
    for (k in seq_along(md)) {
        nm <- names(md)[k] %||% sprintf("item-%d", k)
        if (!nzchar(nm))
            nm <- sprintf("item-%d", k)
        ser <- .serializeMetadataValue(md[[k]], path = path, name = nm,
                                       resources = resources, ...)
        if (!isTRUE(ser$skip))
            annotations[[nm]] <- ser$value
        resources <- ser$resources
    }
    if (length(annotations) == 0L)
        annotations <- NULL
    list(annotations = annotations, resources = resources)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Array-like objects
###

#' @export
#' @importFrom DelayedArray defaultAutoGrid getAutoBPPARAM type
#' @importFrom DuckDBArray writeCoordArray
#' @importFrom S4Vectors head tail
#' @importFrom SparseArray COO_SparseArray
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
         dimension = c("crossed", "unbound"),
         layout = "coord_array",
         ...)
{
    dimension <- match.arg(dimension)

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
                               list(name = indexcols[i], type = "integer")
                           } else {
                               list(name = indexcols[i], type = "integer",
                                    categories =
                                        lapply(seq_along(dimnames_x[[i]]),
                                               function(j) {
                                                   list(value = j,
                                                        label = I(dimnames_x[[i]][[j]]))
                                               }),
                                    categoriesOrdered = TRUE)
                           }
                       }),
                list(.buildFieldSpec(name = datacol, x = vector(type(x)))))

    # Generate foreign key metadata
    if (is.null(indexrefs)) {
        indexrefs <- lapply(seq_along(indexcols), function(i) {
            list(fields = indexcols[i], reference = list(fields = ""))
        })
    }

    schema <- list(fields = fields, foreignKeys = indexrefs)
    resources <- list(list(name = name,
                           path = basename(path),
                           dimension = dimension,
                           layout = layout,
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
            if (inherits(col, "sfc") || is(col, "DuckDBColumn"))
                return(nm)
            if (is.list(col) && length(col) > 0L &&
                all(vapply(col, function(r) is.null(r) || is.raw(r), NA)))
                return(nm)
        }
    }
    NULL
}

.spatialCoordsSchema <- function(x, layout) {
    if (!layout %in% c("spatial_points", "spatial_shapes"))
        return(NULL)
    nms <- colnames(x)
    out <- list()
    for (role in c("x", "y", "z", "instance_id", "geometry")) {
        if (role %in% nms)
            out[[role]] <- role
    }
    geom <- .geometryCol(x)
    if (!is.null(geom))
        out[["geometry"]] <- geom
    if (length(out))
        out
    else
        NULL
}

.buildTableResourceMetadata <-
function(x, path, name, dimension, layout, indexcol, keycol, dimtbl,
         refs = NULL)
{
    has_key <- !is.null(keycol) && keycol %in% colnames(x)

    schema <- list(fields = lapply(colnames(x), function(j) {
        .buildFieldSpec(name = j, x = x[[j]])
    }))

    if (!is.null(indexcol)) {
        schema[["sortOrder"]] <- list(
            list(field = indexcol, direction = "ascending")
        )
    }

    if (has_key && keycol %in% colnames(x)) {
        schema[["primaryKey"]] <- keycol
    } else if (!is.null(indexcol)) {
        schema[["primaryKey"]] <- indexcol
    }

    if (!is.null(refs))
        schema[["foreignKeys"]] <- refs

    if (!is.null(dimtbl) && ncol(dimtbl) > 0L)
        schema[["partitioning"]] <- colnames(dimtbl)

    spatial_coords <- .spatialCoordsSchema(x, layout)
    if (!is.null(spatial_coords))
        schema[["spatialCoords"]] <- spatial_coords

    list(list(name = name,
              path = basename(path),
              dimension = dimension,
              layout = layout,
              format = "parquet",
              mediatype = "application/vnd.apache.parquet",
              schema = schema))
}

.isSubsequentFlatPart <- function(append, part) {
    isTRUE(append) && !is.null(part) && as.integer(part) > 0L
}

.hasLazyGeometryColumn <- function(x) {
    ct <- tryCatch(DuckDBDataFrame::coltypes(x), error = function(e) character())
    any(ct == "geometry") || "geometry" %in% colnames(x)
}

#' @importFrom DuckDBDataFrame dbconn writeDuckDBTableParquet
.writeDuckDBParquet <-
function(x, path, indexcol, keycol, dimtbl, name, dimension, layout,
         refs = NULL, append = FALSE, offset = 0L, part = NULL,
         part_digits = 0L, ...)
{
    geom_write <- identical(layout, "spatial_shapes") || .hasLazyGeometryColumn(x)
    if (geom_write) {
        if (isTRUE(append) || (!is.null(part) && as.integer(part) > 0L))
            stop("flat append ('append', 'part') is not supported for geometry tables")
        if (requireNamespace("DuckDBSpatial", quietly = TRUE)) {
            DuckDBSpatial::enableGeoParquetConversion(dbconn(x))
        }
    }

    result <- writeDuckDBTableParquet(
        x, path = path, indexcol = indexcol, keycol = keycol, dimtbl = dimtbl,
        append = append, offset = offset, part = part, part_digits = part_digits,
        ...)
    if (isTRUE(result$subsequent_part))
        return(invisible(NULL))
    resources <- .buildTableResourceMetadata(
        result$sample_df, path = result$dir, name = name, dimension = dimension,
        layout = layout, indexcol = indexcol, keycol = keycol, dimtbl = dimtbl,
        refs = refs)
    invisible(resources)
}

#' @importFrom arrow Array write_dataset write_parquet
#' @importFrom DuckDBDataFrame arrowType reconcileParquetSchema setupFlatParquetWrite clusterSort
#' @importFrom S4Vectors I
#' @importFrom stats setNames
.writeDataFrameParquet <-
function(x, path, indexcol, keycol, dimtbl, name, dimension, layout,
         refs = NULL, append = FALSE, offset = 0L, part = NULL,
         part_digits = 0L, cluster_by = NULL, ...)
{
    prep <- setupFlatParquetWrite(
        path, append = append, offset = offset, part = part,
        part_digits = part_digits, indexcol = indexcol,
        reconcile_columns = NULL, create = TRUE)
    part <- prep$part
    offset <- prep$offset
    flat_part <- prep$flat_part

    # Cluster the in-memory rows (no SQL ORDER BY on this Arrow path); no-op if NULL. The
    # __index__/__name__ assigned below follow the clustered order.
    if (!is.null(cluster_by))
        x <- clusterSort(x, cluster_by)

    if (is.null(indexcol)) {
        index <- NULL
    } else {
        index <- setNames(list(offset + seq_len(nrow(x))), indexcol)
    }

    rnms <- attr(x, "row.names")
    if (is.null(keycol) || is.integer(rnms)) {
        key <- NULL
    } else {
        key <- setNames(list(rnms), keycol)
    }

    is_sf <- inherits(x, "sf")
    if (is_sf) {
        if (isTRUE(append) || flat_part) {
            stop("flat append ('append', 'part') is not supported for sf objects")
        }
        if (!requireNamespace("DuckDBSpatial", quietly = TRUE)) {
            stop("DuckDBSpatial package required for GeoParquet support; ",
                 "install with BiocManager::install('DuckDBSpatial')")
        }
        x <- do.call(cbind.data.frame, c(index, key, dimtbl, x))
        colnames(x) <- make.unique(colnames(x), sep = "_")
        x <- sf::st_sf(x)
        DuckDBSpatial::writeGeoParquet(x, prep$pq_path,
                         geom = attr(x, "sf_column"), ...)
    } else {
        for (j in seq_along(x)) {
            if (is.list(x[[j]])) {
                x[[j]] <- I(x[[j]])
            }
        }
        x <- do.call(cbind.data.frame, c(index, key, dimtbl, x))
        colnames(x) <- make.unique(colnames(x), sep = "_")

        if (!flat_part && !isTRUE(append)) {
            for (j in seq_along(x)) {
                if (is.integer(x[[j]]) && length(x[[j]]) > 0L) {
                    x[[j]] <- Array$create(x[[j]], type = arrowType(x[[j]]))
                }
            }
        }

        if (isTRUE(append)) {
            int_cols <- colnames(x)[vapply(x, is.integer, logical(1L))]
            if (length(int_cols) > 0L) {
                arrowtypes <- setNames(rep(list(NULL), length(int_cols)),
                                       int_cols)
                reconcileParquetSchema(path, int_cols, arrowtypes)
            }
        }

        if (flat_part) {
            write_parquet(x, prep$pq_path, compression = "zstd",
                          compression_level = 3L, ...)
        } else {
            write_dataset(x, path, format = "parquet", compression = "zstd",
                          compression_level = 3L, ...)
        }
    }

    .buildTableResourceMetadata(
        x, path = path, name = name, dimension = dimension, layout = layout,
        indexcol = indexcol, keycol = keycol, dimtbl = dimtbl, refs = refs)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### data.frame objects
###

#' @export
#' @rdname writeParquet
setMethod("writeParquet", "data.frame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         refs = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("unbound", "sample", "feature", "crossed"),
         layout = "data_frame",
         ...)
{
    dimension <- match.arg(dimension)
    resources <- .writeDataFrameParquet(x, path = path, indexcol = indexcol,
                                        keycol = keycol, dimtbl = dimtbl,
                                        refs = refs, name = name,
                                        dimension = dimension, layout = layout,
                                        append = append, offset = offset,
                                        part = part, part_digits = part_digits,
                                        ...)
    if (.isSubsequentFlatPart(append, part)) {
        return(invisible(NULL))
    }
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
         refs = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("unbound", "sample", "feature", "crossed"),
         layout = "data_frame",
         ...)
{
    dimension <- match.arg(dimension)
    df <- as.data.frame(x, optional = TRUE)

    geom <- .geometryCol(df)
    if (!is.null(geom)) {
        df <- sf::st_as_sf(df, sf_column_name = geom)
    }

    resources <- .writeDataFrameParquet(df, path = path, indexcol = indexcol,
                                        keycol = keycol, dimtbl = dimtbl,
                                        refs = refs, name = name,
                                        dimension = dimension, layout = layout,
                                        append = append, offset = offset,
                                        part = part, part_digits = part_digits,
                                        ...)
    if (.isSubsequentFlatPart(append, part)) {
        return(invisible(NULL))
    }
    invisible(resources)
})

### DuckDBTable objects (lazy SQL COPY TO)
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBTable
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBTable",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("unbound", "sample", "feature", "crossed"),
         layout = "data_frame",
         ...)
{
    dimension <- match.arg(dimension)
    .writeDuckDBParquet(x, path = path, indexcol = indexcol, keycol = keycol,
                        dimtbl = dimtbl, name = name, dimension = dimension,
                        layout = layout, append = append, offset = offset,
                        part = part, part_digits = part_digits, ...)
})

### DuckDBDataFrame objects (lazy write without materialization)
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBDataFrame
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBDataFrame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("unbound", "sample", "feature", "crossed"),
         layout = "data_frame",
         ...)
{
    dimension <- match.arg(dimension)
    .writeDuckDBParquet(x, path = path, indexcol = indexcol, keycol = keycol,
                        dimtbl = dimtbl, name = name, dimension = dimension,
                        layout = layout, append = append, offset = offset,
                        part = part, part_digits = part_digits, ...)
})

### DuckDBTransposedDataFrame objects
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBTransposedDataFrame
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBTransposedDataFrame",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = "crossed",
         layout = "transposed_data_frame",
         ...)
{
    resources <- callGeneric(t(x), path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             append = append, offset = offset, part = part,
                             part_digits = part_digits, name = name,
                             dimension = dimension, layout = layout, ...)
    invisible(resources)
})

### DuckDBSelfHits objects
###

#' @export
#' @importClassesFrom DuckDBDataFrame DuckDBSelfHits
#' @importFrom S4Vectors nnode
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBSelfHits",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("unbound", "sample", "feature"),
         layout = "graph_edges",
         ...)
{
    dimension <- match.arg(dimension)
    resources <- callGeneric(x@frame, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             append = append, offset = offset, part = part,
                             part_digits = part_digits, name = name,
                             dimension = dimension, layout = layout, ...)
    if (is.null(resources))
        return(invisible(NULL))
    schema <- resources[[length(resources)]][["schema"]]
    schema <- .addGraphMetadata(schema, nnode(x))
    resources[[length(resources)]][["schema"]] <- schema
    invisible(resources)
})

### DuckDBGRanges objects
###

#' @export
#' @importClassesFrom DuckDBGRanges DuckDBGRanges
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBGRanges",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("feature", "unbound"),
         layout = "genomic_ranges",
         ...)
{
    dimension <- match.arg(dimension)
    resources <- callGeneric(x@frame, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             append = append, offset = offset, part = part,
                             part_digits = part_digits, name = name,
                             dimension = dimension, layout = layout, ...)
    if (is.null(resources))
        return(invisible(NULL))
    schema <- resources[[length(resources)]][["schema"]]
    schema <- .addGenomicMetadata(schema)
    resources[[length(resources)]][["schema"]] <- schema
    invisible(resources)
})

### DuckDBGRangesList objects
###

#' @export
#' @importClassesFrom DuckDBGRanges DuckDBGRangesList
#' @rdname writeParquet
setMethod("writeParquet", "DuckDBGRangesList",
function(x,
         path,
         indexcol = "__index__",
         keycol = "__name__",
         dimtbl = NULL,
         append = FALSE,
         offset = 0L,
         part = NULL,
         part_digits = 0L,
         name = basename(path),
         dimension = c("feature", "unbound"),
         layout = "genomic_ranges_list",
         ...)
{
    dimension <- match.arg(dimension)
    resources <- callGeneric(x@frame, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             append = append, offset = offset, part = part,
                             part_digits = part_digits, name = name,
                             dimension = dimension, layout = layout, ...)
    if (is.null(resources))
        return(invisible(NULL))
    schema <- resources[[length(resources)]][["schema"]]
    schema <- .addGenomicMetadata(schema)
    resources[[length(resources)]][["schema"]] <- schema
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
         dimension = "crossed",
         layout = "transposed_data_frame",
         ...)
{
    resources <- callGeneric(t(x), path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, dimension = dimension,
                             layout = layout, ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list helper functions
###

.writeListParquet <- function(x, path, dimension, ...)
{
    layout <- list(...)[["layout"]]
    if (is.null(layout)) {
        prefix <- ""
    } else {
        prefix <- switch(layout,
                         coord_array =,
                         data_frame =,
                         transposed_data_frame = "assay_",
                         nested_data_frame     = "table_",
                         graph_edges           = "graph_",
                         spatial_points        = "points_",
                         spatial_shapes        = "shapes_",
                         stop("unsupported layout: ", layout))
        if (dimension != "crossed") {
            prefix <- sprintf("%s_%s", dimension, prefix)
        }
    }

    resources <- list()
    for (i in seq_along(x)) {
        x_i <- getListElement(x, i)
        nms_i <- names(x)[i] %||% sprintf("item-%d", i)
        path_i <- sprintf("%s%s", prefix, nms_i)
        resources <- c(resources,
                       writeParquet(x_i, path = file.path(path, path_i),
                                    name = nms_i, dimension = dimension, ...))
    }
    resources
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### list objects
###

#' @export
#' @importFrom S4Vectors getListElement
#' @rdname writeParquet
setMethod("writeParquet", "list",
function(x,
         path,
         dimension = c("unbound", "sample", "feature", "crossed"),
         ...)
{
    dimension <- match.arg(dimension)
    resources <- .writeListParquet(x, path = path, dimension = dimension, ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List objects
###

#' @export
#' @importFrom S4Vectors getListElement
#' @rdname writeParquet
setMethod("writeParquet", "List",
function(x,
         path,
         dimension = c("unbound", "sample", "feature", "crossed"),
         ...)
{
    dimension <- match.arg(dimension)
    resources <- .writeListParquet(x, path = path, dimension = dimension, ...)
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
         dimension = c("feature", "unbound"),
         layout = "genomic_ranges",
         ...)
{
    dimension <- match.arg(dimension)

    # Write Data Table
    df <- as.data.frame(x, optional = TRUE)
    df[["width"]] <- NULL
    resources <- callGeneric(df, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, dimension = dimension,
                             layout = layout, ...)

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
         dimension = c("feature", "unbound"),
         layout = "genomic_ranges_list",
         ...)
{
    dimension <- match.arg(dimension)

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
    resources <- callGeneric(df, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, dimension = dimension,
                             layout = layout, ...)

    # Add domain-specific metadata to schema for ranges
    schema <- resources[[length(resources)]][["schema"]]

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
         dimension = c("unbound", "sample", "feature"),
         layout = "graph_edges",
         ...)
{
    dimension <- match.arg(dimension)
    df <- as.data.frame(x, optional = TRUE)
    resources <- callGeneric(df, path = path, indexcol = indexcol,
                             keycol = keycol, dimtbl = dimtbl,
                             name = name, dimension = dimension,
                             layout = layout, ...)

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
         dimension = c("sample", "feature"),
         ...)
{
    dimension <- match.arg(dimension)
    x <- lapply(x, function(y) if (is(y, "DualSubset")) y@hits  else y)
    resources <- .writeListParquet(x, path = path, dimension = dimension,
                                   layout = "graph_edges", ...)
    invisible(resources)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assays objects
###

#' @export
#' @importClassesFrom S4Vectors DataFrame TransposedDataFrame
#' @importClassesFrom SummarizedExperiment Assays
#' @importFrom S4Vectors getListElement
#' @rdname writeParquet
setMethod("writeParquet", "Assays",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         indexrefs = NULL,
         datacol = "value",
         grid = defaultAutoGrid(COO_SparseArray(dim(x[[1L]]))),
         grid_suffix = "group__",
         dimension = "crossed",
         ...)
{
    # Write transpose of the assays
    indexcols <- rev(indexcols)
    indexrefs <- rev(indexrefs)
    grid <- t(grid)
    resources <- list()
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- t(getListElement(x, i))
        path_i <- paste0("assay_", nms_i)
        if (is(x_i, "DataFrame")) {
            resources <- c(resources,
                           callGeneric(x_i, path = file.path(path, path_i),
                                       indexcol = indexcols[1L], name = nms_i,
                                       dimension = dimension, ...))
        } else if (is(x_i, "TransposedDataFrame")) {
            resources <- c(resources,
                           callGeneric(x_i, path = file.path(path, path_i),
                                       indexcol = indexcols[2L], name = nms_i,
                                       dimension = dimension, ...))
        } else {
            # Array-like object
            dimnames(x_i) <- NULL
            resources <- c(resources,
                           callGeneric(x_i, path = file.path(path, path_i),
                                       indexcols = indexcols,
                                       indexrefs = indexrefs, datacol = datacol,
                                       grid = grid, grid_suffix = grid_suffix,
                                       name = nms_i, dimension = dimension,
                                       ...))
        }
    }

    invisible(resources)
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
         package = list(model = ifelse(is(x, "RangedSummarizedExperiment"),
                                       "ranged_summarized_experiment",
                                       "summarized_experiment"),
                        resources = list()),
         ...)
{
    # Make dimnames unique
    rownames(x) <- make.unique(rownames(x) %||% as.character(seq_len(nrow(x))),
                               sep = "_")
    colnames(x) <- make.unique(colnames(x) %||% as.character(seq_len(ncol(x))),
                               sep = "_")

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
                    indexcol = indexcols[1L], dimtbl = dimtbls[[1L]],
                    dimension = "feature", ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Sample Data
    resources <-
        callGeneric(colData(x), path = file.path(path, "samples"),
                    indexcol = indexcols[2L], dimtbl = dimtbls[[2L]],
                    dimension = "sample", ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Assay Data
    indexrefs <- list(list(fields = indexcols[1L],
                           reference = list(fields = indexcols[1L],
                                            resource = "features")),
                      list(fields = indexcols[2L],
                           reference = list(fields = indexcols[2L],
                                            resource = "samples")))
    resources <- callGeneric(x@assays, path = path,
                             indexcols = indexcols, indexrefs = indexrefs,
                             grid = grid, grid_suffix = grid_suffix,
                             dimension = "crossed", ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Metadata — recursive JSON vs Parquet dispatch
    ser <- .serializeMetadata(x, path = path, ...)
    package[["resources"]] <- c(package[["resources"]], ser$resources)
    package[["annotations"]] <- ser[["annotations"]]

    # Declare the Frictionless profile
    package <- c(list("$schema" = .BIOCDUCKDB_PROFILE), package)

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Additional Dimension Tables
###

#' @importFrom jsonlite write_json
#' @importFrom S4Vectors getListElement
.writeParquetDimTables <-
function(x,
         path,
         dimension = c("sample", "feature"),
         indexcol = "__index__",
         ...)
{
    dimension <- match.arg(dimension)
    for (i in seq_along(x)) {
        x_i <- getListElement(x, i)
        if (!is(x_i, "DataFrame"))
            stop("Unsupported object type: ", class(x_i))
    }
    resources <- .writeListParquet(x, path = path, dimension = dimension,
                                   layout = "nested_data_frame",
                                   indexcol = indexcol, ...)
    invisible(resources)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SingleCellExperiment objects
###

#' @export
#' @importClassesFrom S4Vectors DFrame
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment ExperimentList
#' @importFrom S4Vectors I
#' @importFrom SingleCellExperiment altExps colPairs mainExpName
#' @importFrom SingleCellExperiment reducedDims reducedDimNames rowPairs
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
#' @rdname writeParquet
setMethod("writeParquet", "SingleCellExperiment",
function(x,
         path,
         indexcols = c("__feature__", "__sample__"),
         package = list(model = "single_cell_experiment",
                        resources = list()),
         ...)
{
    # Make dimnames unique
    rownames(x) <- make.unique(rownames(x) %||% as.character(seq_len(nrow(x))),
                               sep = "_")
    colnames(x) <- make.unique(colnames(x) %||% as.character(seq_len(ncol(x))),
                               sep = "_")

    # Unnest rowData
    rdata <- rowData(x)
    if (is(rdata, "DFrame")) {
        nested <- unlist(lapply(rdata, is, "DFrame"))
        nested <- names(nested)[nested]
        for (j in nested) {
            rowTable(x, j) <- rdata[[j]]
            rdata[[j]] <- NULL
        }
        if (length(nested)) {
            rowData(x) <- rdata
        }
    }

    # Unnest colData
    cdata <- colData(x)
    if (is(cdata, "DFrame")) {
        nested <- unlist(lapply(cdata, is, "DFrame"))
        nested <- names(nested)[nested]
        for (j in nested) {
            colTable(x, j) <- cdata[[j]]
            cdata[[j]] <- NULL
        }
        if (length(nested)) {
            colData(x) <- cdata
        }
    }

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
                                 indexcol = indexcols[1L],
                                 dimension = "feature",
                                 layout = "embedding_table",
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
                                 indexcol = indexcols[2L],
                                 dimension = "sample",
                                 layout = "embedding_table",
                                 ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Alternative Experiments
    exps <- ExperimentList(altExps(x))
    if (length(exps)) {
        resources <- callGeneric(exps, path = path, indexcols = indexcols, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Row Tables
    rtables <- rowTables(x)
    if (length(rtables)) {
        resources <- .writeParquetDimTables(rtables, path = path,
                                            indexcol = indexcols[1L],
                                            dimension = "feature", ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Column Tables
    ctables <- colTables(x)
    if (length(ctables)) {
        resources <- .writeParquetDimTables(ctables, path = path,
                                            indexcol = indexcols[2L],
                                            dimension = "sample", ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Row Pairs
    rpairs <- rowPairs(x, asSparse = FALSE)
    if (length(rpairs)) {
        resources <- .writeParquetGraphs(rpairs, path = path,
                                         dimension = "feature", ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Column Pairs
    cpairs <- colPairs(x, asSparse = FALSE)
    if (length(cpairs)) {
        resources <- .writeParquetGraphs(cpairs, path = path,
                                         dimension = "sample", ...)
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
         ...)
{
    resources <- list()
    for (i in seq_along(x)) {
        nms_i <- names(x)[i]
        x_i <- x[[i]]
        path_i <- paste0("experiment_", nms_i)
        if (is(x_i, "SummarizedExperiment")) {
            res_i <- list(list(name = nms_i, path = path_i,
                               dimension = "crossed",
                               layout = "nested_experiment"))
           callGeneric(x_i, path = file.path(path, path_i),
                       indexcols = indexcols, ...)
        } else {
            # Array-like object
            res_i <-
                callGeneric(x_i, path = file.path(path, path_i),
                            indexcols = indexcols, datacol = "value",
                            grid_suffix = "group__", name = nms_i,
                            dimension = "crossed", layout = "coord_array", ...)
        }
        resources <- c(resources, res_i)
    }
    invisible(resources)
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
         package = list(model = "multi_assay_experiment",
                        resources = list()),
         ...)
{
    # Subject Data
    subj_res <- callGeneric(colData(x), path = file.path(path, "subjects"),
                            dimension = "sample", ...)
    package[["resources"]] <- c(package[["resources"]], subj_res)

    # Sample Map
    sm <- sampleMap(x)
    refs <- NULL
    subj <- Find(function(r) identical(r[["name"]], "subjects"), subj_res)
    if (!is.null(subj) && !is.null(rownames(colData(x))) &&
        "primary" %in% colnames(sm)) {
        refs <- list(list(fields = "primary",
                          reference = list(
                              fields = subj[["schema"]][["primaryKey"]],
                              resource = "subjects")))
    }
    resources <- callGeneric(sm, path = file.path(path, "sample_map"),
                             dimension = "unbound", refs = refs, ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    # Experiment Data
    resources <- callGeneric(experiments(x), path = path, indexcols = indexcols,
                             ...)
    package[["resources"]] <- c(package[["resources"]], resources)

    ser <- .serializeMetadata(x, path = path, ...)
    package[["resources"]] <- c(package[["resources"]], ser$resources)
    package[["annotations"]] <- ser[["annotations"]]

    # Declare the Frictionless profile
    package <- c(list("$schema" = .BIOCDUCKDB_PROFILE), package)

    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(NULL)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PointsLayerList objects
###

#' @export
#' @importClassesFrom MultiAssaySpatialExperiment PointsLayerList
#' @rdname writeParquet
setMethod("writeParquet", "PointsLayerList", function(x, path, ...)
{
    resources <- .writeListParquet(x, path = path, dimension = "sample",
                                   layout = "spatial_points", ...)
    invisible(resources)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ShapesLayerList objects
###

#' @export
#' @importClassesFrom MultiAssaySpatialExperiment ShapesLayerList
#' @rdname writeParquet
setMethod("writeParquet", "ShapesLayerList", function(x, path, ...)
{
    resources <- .writeListParquet(x, path = path, dimension = "sample",
                                   layout = "spatial_shapes", ...)
    invisible(resources)
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
         package = list(model = "multi_assay_spatial_experiment",
                        resources = list()),
         ...)
{
    # Points
    pts <- MultiAssaySpatialExperiment::spatialPoints(x)
    if (length(pts) > 0L) {
        resources <- callGeneric(pts, path = path, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Shapes
    shps <- MultiAssaySpatialExperiment::spatialShapes(x)
    if (length(shps) > 0L) {
        resources <- callGeneric(shps, path = path, ...)
        package[["resources"]] <- c(package[["resources"]], resources)
    }

    # Images path references
    imgs <- MultiAssaySpatialExperiment::spatialImages(x)
    if (length(imgs) > 0L) {
        imgs_dir <- file.path(path, "images")
        dir.create(imgs_dir, showWarnings = FALSE, recursive = TRUE)
        img_resources <- list()
        for (i in seq_along(imgs)) {
            nm <- names(imgs)[i]
            el_path <- file.path(imgs_dir, nm)
            .writeRasterRef(imgs[[i]], el_path, nm, type = "image")
            img_resources <- c(img_resources, list(list(
                name = paste0("sample_images_", nm),
                path = file.path("images", nm),
                dimension = "sample",
                layout = "spatial_raster_ref",
                format = "json",
                mediatype = "application/json")))
        }
        package[["resources"]] <- c(package[["resources"]], img_resources,
                                    list(list(name = "images",
                                              path = "images")))
    }

    # Labels path references
    lbls <- MultiAssaySpatialExperiment::spatialLabels(x)
    if (length(lbls) > 0L) {
        lbls_dir <- file.path(path, "labels")
        dir.create(lbls_dir, showWarnings = FALSE, recursive = TRUE)
        lab_resources <- list()
        for (i in seq_along(lbls)) {
            nm <- names(lbls)[i]
            el <- lbls[[i]]
            el_path <- file.path(lbls_dir, nm)
            if (is.matrix(el) || is.array(el)) {
                coo_path <- file.path(el_path, "coord")
                .matrixToCoordLabel(el, coo_path)
                lab_resources <- c(lab_resources, list(list(
                    name = paste0("sample_labels_", nm),
                    path = file.path("labels", nm, "coord"),
                    dimension = "sample",
                    layout = "spatial_label_coord",
                    format = "parquet",
                    mediatype = "application/vnd.apache.parquet")))
            } else {
                .writeRasterRef(el, el_path, nm, type = "label")
                lab_resources <- c(lab_resources, list(list(
                    name = paste0("sample_labels_", nm),
                    path = file.path("labels", nm),
                    dimension = "sample",
                    layout = "spatial_raster_ref",
                    format = "json",
                    mediatype = "application/json")))
            }
        }
        package[["resources"]] <- c(package[["resources"]], lab_resources,
                                    list(list(name = "labels",
                                              path = "labels")))
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

            resources_img <- callGeneric(img_data_copy,
                                         path = file.path(path, "img_data"),
                                         dimension = "unbound", ...)
        } else {
            resources_img <- callGeneric(img_data,
                                         path = file.path(path, "img_data"),
                                         dimension = "unbound", ...)
        }
        package[["resources"]] <- c(package[["resources"]], resources_img)
    }

    # Spatial Map
    spatial_map <- MultiAssaySpatialExperiment::spatialMap(x)
    if (!is.null(spatial_map) && nrow(spatial_map) > 0L) {
        resources_sm <- callGeneric(spatial_map, path = file.path(path, "spatial_map"),
                                    dimension = "unbound", ...)
        package[["resources"]] <- c(package[["resources"]], resources_sm)
    }

    callNextMethod(x, path = path, indexcols = indexcols, package = package, ...)

    invisible(NULL)
})
