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
#' \code{datapackage.json} file that describes the structure and class of the
#' Bioconductor object following the Frictionless Data Package specification
#' and the written objects must have been created with
#' \code{\link{writeParquet}}.
#' @param metadata A list containing metadata about the object structure. If
#' \code{NULL} (default), the metadata is read from
#' \code{file.path(path, "datapackage.json")}. The metadata should contain
#' information about the object type, resources, schemas, and any additional
#' components specific to the object type following the Frictionless Data
#' Package specification and the written objects must have been created with
#' \code{\link{writeParquet}}.
#' @param class A character string specifying the class of Bioconductor object to
#' read. If \code{NULL} (default), the class is determined from the metadata.
#' Supported types include:
#' \itemize{
#'   \item \code{"assays"} - \code{Assays} objects
#'   \item \code{"genomic_ranges"} - \code{GenomicRanges} objects
#'   \item \code{"genomic_ranges_list"} - \code{GenomicRangesList} objects
#'   \item \code{"graph_edges"} - \code{DuckDBSelfHits} objects
#'   \item \code{"summarized_experiment"} - Basic \code{SummarizedExperiment} objects
#'   \item \code{"ranged_summarized_experiment"} - \code{RangedSummarizedExperiment} objects
#'   \item \code{"single_cell_experiment"} - \code{SingleCellExperiment} objects
#'   \item \code{"experiment_list"} - \code{ExperimentList} objects
#'   \item \code{"multi_assay_experiment"} - \code{MultiAssayExperiment} objects
#'   \item \code{"multi_assay_spatial_experiment"} - \code{MultiAssaySpatialExperiment} objects
#' }
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
#' functionality by also reading reduced dimensions, alternative experiments, and
#' column/row pairings from their respective parquet files. Pairwise graphs stored
#' in \code{sample_graphs/} and \code{feature_graphs/} are reconstructed as
#' \linkS4class{DuckDBDualSubset} objects wrapping \linkS4class{DuckDBSelfHits}.
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
#'   \item{\code{Assays}}{
#'     \code{Assays} objects are read from the \code{assays/} subdirectory.
#'   }
#'   \item{\code{GenomicRanges}}{
#'     \code{GenomicRanges} objects are read with genomic coordinates (seqnames,
#'     start, end, strand) and metadata columns.
#'   }
#'   \item{\code{GenomicRangesList}}{
#'     \code{GenomicRangesList} objects are read with genomic coordinates
#'     (seqnames, start, end, strand) and metadata columns.
#'   }
#'   \item{\code{SummarizedExperiment}}{
#'     Basic multi-assay genomic experiments with feature and sample metadata
#'     from the \code{features/} and \code{samples/} subdirectories.
#'   }
#'   \item{\code{RangedSummarizedExperiment}}{
#'     \code{SummarizedExperiment} with genomic ranges for features from the
#'     \code{features/} subdirectory.
#'   }
#'   \item{\code{SingleCellExperiment}}{
#'     Single-cell genomic experiments with reduced dimensions, alternative
#'     experiments, and pairwise graphs (\code{colPairs}/\code{rowPairs}) from the
#'     \code{embeddings/}, \code{modalities/}, \code{sample_graphs/}, and
#'     \code{feature_graphs/} subdirectories.
#'   }
#'   \item{\code{DuckDBSelfHits}}{
#'     Graph edge lists with \code{from}, \code{to} columns and optional metadata.
#'     Node count (\code{nnode}) and column names are stored in the schema
#'     \code{graphCoords} property and used to reconstruct the object.
#'   }
#'   \item{\code{ExperimentList}}{
#'     \code{ExperimentList} objects are read from the \code{experiments/} subdirectory.
#'   }
#'   \item{\code{MultiAssayExperiment}}{
#'     Multi-experiment studies that link experiments from the \code{experiments/}
#'     subdirectory.
#'   }
#'   \item{\code{MultiAssaySpatialExperiment}}{
#'     Multi-assay spatial experiments with points, shapes, imgData, and spatialMap
#'     from the \code{points/}, \code{shapes/}, \code{img_data/}, \code{spatial_map/}
#'     subdirectories. Requires the \code{MultiAssaySpatialExperiment} package.
#'   }
#' }
#'
#' @section Frictionless Data Package Metadata:
#' This function reads data packages that follow the Frictionless Data Package
#' specification, parsing the \code{datapackage.json} metadata file to understand
#' the data structure. The metadata provides:
#' \itemize{
#'   \item Resource definitions with paths and types
#'   \item Schema information including field types and constraints
#'   \item Primary key definitions for data integrity
#'   \item Semantic annotations for genomic data (sequence names, coordinates)
#'   \item Object type information for proper reconstruction
#' }
#' This standardized approach ensures consistent data interpretation and
#' facilitates interoperability with other tools that support the Frictionless
#' Data Package specification.
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
#' }
#'
#' @include fieldtypes.R
#'
#' @export
#' @importFrom jsonlite read_json
#' @rdname readParquet
readParquet <-
function(path,
         metadata = read_json(file.path(path, "datapackage.json"),
                              simplifyVector = TRUE,
                              simplifyDataFrame = FALSE,
                              simplifyMatrix = FALSE),
         class = metadata[["class"]],
         ...)
{
    switch(class,
           "data_frame"              = .readParquetDataFrame(path, metadata, ...),
           "transposed_data_frame"   = .readParquetTransposedDataFrame(path, metadata, ...),
           "array"                   = .readParquetArray(path, metadata, ...),
           "data_package"            = .readParquetDataPackage(path, metadata, ...),
           "genomic_ranges"          = .readParquetGenomicRanges(path, metadata, ...),
           "genomic_ranges_list"     = .readParquetGenomicRangesList(path, metadata, ...),
           "graph_edges"             = .readParquetGraphEdges(path, metadata, ...),
           "assays"                  = .readParquetAssays(path, metadata, ...),
           "ranged_summarized_experiment" =,
           "summarized_experiment"   =,
           "single_cell_experiment"  = .readParquetSE(path, metadata, ...),
           "experiment_list"         = .readParquetExps(path, metadata, ...),
           "multi_assay_experiment"  = .readParquetMAE(path, metadata, ...),
           "multi_assay_spatial_experiment" = .readParquetMASE(path, metadata, ...),
           stop("unsupported Bioconductor parquet layout"))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

.schema_keycols <- function(schema) {
    if (is.null(schema[["foreignKeys"]])) {
        schema[["sortOrder"]][[1L]][["field"]] %||% schema[["primaryKey"]]
    } else {
        unlist(lapply(schema[["foreignKeys"]], `[[`, "fields"))
    }
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

.schema_genomic <- function(schema, role) {
    schema[["genomicCoords"]][[role]]
}

.schema_graph <- function(schema, role) {
    schema[["graphEdges"]][[role]]
}

.schema_partitions <- function(schema) {
    schema[["partitioning"]]
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
    DuckDBDataFrame(path, datacols = datacols, keycol = keycol)
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
### Generic Data Packages
###

.readParquetDataPackage <- function(path, package, ...) {
    lst <- lapply(package[["resources"]], function(x) {
        readParquet(file.path(path, x[["path"]]), metadata = x, ...)
    })
    names(lst) <- sapply(package[["resources"]], `[[`, "name")
    lst
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assays objects
###

.readParquetAssays <- function(path, package, keycols, dimtbls, ...) {
    # Assays are written transposed, so we need to transpose them back
    assays <- lapply(package[["resources"]], function(x) {
        t(readParquet(file.path(path, x[["path"]]), metadata = x,
                      keycols = rev(keycols), dimtbls = dimtbls, ...))
    })
    names(assays) <- sapply(package[["resources"]], `[[`, "name")
    assays
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SummarizedExperiment objects
###

#' @importFrom DuckDBDataFrame DuckDBDataFrame
#' @importFrom jsonlite read_json
#' @importFrom IRanges DataFrameList
#' @importFrom S4Vectors endoapply
#' @importFrom SingleCellExperiment SingleCellExperiment colPair<- rowPair<-
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment SummarizedExperiment
.readParquetSE <- function(path, package, class = package[["class"]], ...) {
    resources <- package[["resources"]]
    names(resources) <- sapply(package[["resources"]], `[[`, "name")

    # Dimension Dictionaries
    dimdicts <- list(features = resources[["features"]][["schema"]],
                     samples = resources[["samples"]][["schema"]])

    # Dimension Tables
    fun <- function(path, i) {
        keycol <- .schema_keycols(dimdicts[[i]])
        pkey <- dimdicts[[i]][["primaryKey"]]
        name <- if (!identical(pkey, keycol)) pkey else NULL
        partitions <- .schema_partitions(dimdicts[[i]])
        df <- DuckDBDataFrame(file.path(path, resources[[i]][["path"]]),
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
    dimtbls[["dimtbls"]] <-
        as(c(fun(path, "features"), fun(path, "samples")), "DataFrameList")
    dimkeycols <- lapply(dimtbls[["dimtbls"]],
                         function(x) setNames(seq_len(nrow(x)), rownames(x)))
    dimtbls[["dimtbls"]] <- endoapply(dimtbls[["dimtbls"]], function(x) {
        rownames(x) <- NULL
        x
    })

    # Feature Data
    schema <- dimdicts[["features"]]
    keycol <- dimkeycols[1L]

    feature_data <- NULL
    feature_ranges <- NULL
    if (is.null(schema[["genomicCoords"]])) {
        fullpath <- file.path(path, resources[["features"]][["path"]])
        feature_data <- readParquet(fullpath,
                                    metadata = resources[["features"]],
                                    keycol = keycol)
    } else {
        fullpath <- file.path(path, resources[["features"]][["path"]])
        feature_ranges <- readParquet(fullpath,
                                      metadata = resources[["features"]],
                                      keycol = keycol)
    }

    # Sample Data
    keycol <- dimkeycols[2L]
    fullpath <- file.path(path, resources[["samples"]][["path"]])
    samples <- readParquet(fullpath, metadata = resources[["samples"]],
                           keycol = keycol)

    # Assays
    fullpath <- file.path(path, resources[["assays"]][["path"]])
    assays <- readParquet(fullpath, keycols = dimkeycols, dimtbls = dimtbls)

    # Metadata
    metadata <- package[["annotations"]]

    # SummarizedExperiment
    if (class %in% c("ranged_summarized_experiment", "summarized_experiment")) {
        se <- SummarizedExperiment(assays,
                                   rowData = feature_data,
                                   rowRanges = feature_ranges,
                                   colData = samples,
                                   metadata = metadata)
    } else if (class == "single_cell_experiment") {
        # Embeddings
        if (is.null(resources[["embeddings"]])) {
            rdims <- list()
        } else {
            keycol <- dimkeycols[2L]
            fullpath <- file.path(path, resources[["embeddings"]][["path"]])
            rdims <- readParquet(fullpath, metadata = resources[["embeddings"]],
                                 keycol = keycol)
            rdims <- as.list(rdims)
        }

        # Modalities
        if (is.null(resources[["modalities"]])) {
            alts <- list()
        } else {
            dirpath <- file.path(path, resources[["modalities"]][["path"]])
            pkg <- read_json(file.path(dirpath, "datapackage.json"),
                             simplifyVector = TRUE,
                             simplifyDataFrame = FALSE,
                             simplifyMatrix = FALSE)
            alts <- lapply(pkg[["resources"]], function(x) {
                readParquet(file.path(dirpath, x[["path"]]), ...)
            })
            names(alts) <- sapply(pkg[["resources"]], `[[`, "name")
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

        # Row Graphs
        if (!is.null(resources[["feature_graphs"]])) {
            dirpath <- file.path(path, resources[["feature_graphs"]][["path"]])
            pkg <- read_json(file.path(dirpath, "datapackage.json"),
                             simplifyVector = TRUE,
                             simplifyDataFrame = FALSE,
                             simplifyMatrix = FALSE)
            for (i in seq_along(pkg[["resources"]])) {
                res <- pkg[["resources"]][[i]]
                hits <- readParquet(file.path(dirpath, res[["path"]]),
                                    metadata = res, ...)
                rowPair(se, res[["name"]]) <- hits
            }
        }

        # Column Graphs
        if (!is.null(resources[["sample_graphs"]])) {
            dirpath <- file.path(path, resources[["sample_graphs"]][["path"]])
            pkg <- read_json(file.path(dirpath, "datapackage.json"),
                             simplifyVector = TRUE,
                             simplifyDataFrame = FALSE,
                             simplifyMatrix = FALSE)
            for (i in seq_along(pkg[["resources"]])) {
                res <- pkg[["resources"]][[i]]
                hits <- readParquet(file.path(dirpath, res[["path"]]),
                                    metadata = res, ...)
                colPair(se, res[["name"]]) <- hits
            }
        }
    } else {
        stop("unsupported Bioconductor parquet layout")
    }

    se
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ExperimentList objects
###

#' @importFrom jsonlite read_json
.readParquetExps <- function(path, package, ...) {
    exps <- lapply(package[["resources"]], function(x) {
                       fullpath <- file.path(path, x[["path"]])
                       jsonfile <- file.path(fullpath, "datapackage.json")
                       if (file.exists(jsonfile)) {
                           x <- read_json(jsonfile,
                                          simplifyVector = TRUE,
                                          simplifyDataFrame = FALSE,
                                          simplifyMatrix = FALSE)
                       }
                       readParquet(fullpath, metadata = x, ...)
                   })
    names(exps) <- sapply(package[["resources"]], `[[`, "name")
    exps
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### MultiAssayExperiment objects
###

#' @importFrom jsonlite read_json
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom stats setNames
.readParquetMAE <- function(path, package, ...) {
    resources <- package[["resources"]]
    names(resources) <- sapply(package[["resources"]], `[[`, "name")

    # Experiments
    fullpath <- file.path(path, resources[["experiments"]][["path"]])
    experiments <- readParquet(fullpath, ...)

    # Subject Data
    sample_schema <- resources[["subjects"]][["schema"]]
    index <- .schema_keycols(sample_schema)
    pkey <- sample_schema[["primaryKey"]]
    name <- if (!is.null(pkey) && !identical(pkey, index)) pkey else NULL
    fullpath <- file.path(path, resources[["subjects"]][["path"]])
    keycol <- readParquet(fullpath, metadata = resources[["subjects"]],
                          datacols = name, keycol = index)
    keycol <- as.list(as.data.frame(keycol, optional = TRUE))
    subjects <- readParquet(fullpath, metadata = resources[["subjects"]],
                            keycol = keycol)

    # Sample Map
    index <- .schema_keycols(resources[["sample_map"]][["schema"]])
    fullpath <- file.path(path, resources[["sample_map"]][["path"]])
    sample_map <- readParquet(fullpath, metadata = resources[["sample_map"]],
                              keycol = index)
    sample_map <- as.data.frame(sample_map, optional = TRUE)
    sample_map[[1L]] <- factor(sample_map[[1L]], levels = names(experiments))

    # Metadata
    metadata <- package[["annotations"]] %||% list()

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

    # Points (optional) - materialize for MASE validation
    points <- MultiAssaySpatialExperiment::PointsLayerList()
    if (!is.null(resources[["points"]])) {
        fullpath <- file.path(path, resources[["points"]][["path"]])
        pkg <- read_json(file.path(fullpath, "datapackage.json"),
                         simplifyVector = TRUE,
                         simplifyDataFrame = FALSE,
                         simplifyMatrix = FALSE)
        pts_list <- lapply(pkg[["resources"]], function(r) {
            x <- readParquet(file.path(fullpath, r[["path"]]), metadata = r, ...)
            as(x, "DataFrame")
        })
        names(pts_list) <- sapply(pkg[["resources"]], `[[`, "name")
        points <- MultiAssaySpatialExperiment::PointsLayerList(pts_list)
    }

    # Shapes (optional) - materialize for MASE validation
    shapes <- MultiAssaySpatialExperiment::ShapesLayerList()
    if (!is.null(resources[["shapes"]])) {
        fullpath <- file.path(path, resources[["shapes"]][["path"]])
        pkg <- read_json(file.path(fullpath, "datapackage.json"),
                         simplifyVector = TRUE,
                         simplifyDataFrame = FALSE,
                         simplifyMatrix = FALSE)
        shp_list <- lapply(pkg[["resources"]], function(r) {
            x <- readParquet(file.path(fullpath, r[["path"]]), metadata = r, ...)
            as(x, "DataFrame")
        })
        names(shp_list) <- sapply(pkg[["resources"]], `[[`, "name")
        shapes <- MultiAssaySpatialExperiment::ShapesLayerList(shp_list)
    }

    # Images (optional) - path references
    images <- MultiAssaySpatialExperiment::RasterLayerList()
    if (!is.null(resources[["images"]])) {
        fullpath <- file.path(path, resources[["images"]][["path"]])
        if (dir.exists(fullpath)) {
            json_files <- list.files(fullpath, pattern = "\\.json$", full.names = TRUE)
            if (length(json_files) > 0L) {
                img_list <- lapply(json_files, function(jf) {
                    meta <- read_json(jf, simplifyVector = TRUE)
                    ## TODO: construct StoredSpatialImage or similar from path ref
                    ## For now, store metadata as-is
                    meta
                })
                names(img_list) <- sub("\\.json$", "", basename(json_files))
                images <- MultiAssaySpatialExperiment::RasterLayerList(img_list)
            }
        }
    }

    # Labels (optional) - path references
    labels <- MultiAssaySpatialExperiment::RasterLayerList()
    if (!is.null(resources[["labels"]])) {
        fullpath <- file.path(path, resources[["labels"]][["path"]])
        if (dir.exists(fullpath)) {
            json_files <- list.files(fullpath, pattern = "\\.json$", full.names = TRUE)
            if (length(json_files) > 0L) {
                lbl_list <- lapply(json_files, function(jf) {
                    meta <- read_json(jf, simplifyVector = TRUE)
                    ## TODO: construct appropriate label/mask object from path ref
                    meta
                })
                names(lbl_list) <- sub("\\.json$", "", basename(json_files))
                labels <- MultiAssaySpatialExperiment::RasterLayerList(lbl_list)
            }
        }
    }

    # imgData (optional) - must materialize
    if (!is.null(resources[["img_data"]])) {
        fullpath <- file.path(path, resources[["img_data"]][["path"]])
        img_data <- readParquet(fullpath, metadata = resources[["img_data"]], ...)
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
        index <- .schema_keycols(resources[["spatial_map"]][["schema"]])
        fullpath <- file.path(path, resources[["spatial_map"]][["path"]])
        x <- readParquet(fullpath, metadata = resources[["spatial_map"]], keycol = index, ...)
        spatial_map <- as.data.frame(x, optional = TRUE)
        spatial_map <- S4Vectors::DataFrame(spatial_map)
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
