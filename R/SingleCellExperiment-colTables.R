#' Sample table methods
#'
#' @description
#' Methods to get or set nested sample-level tables in a
#' \linkS4class{SingleCellExperiment} object. These tables enable storage of
#' multi-column relational data associated with samples (cells/perturbations),
#' which would otherwise require nested DataFrames in \code{colData()}. By
#' extracting nested tables into a separate slot, they can be serialized to
#' independent Parquet tables for efficient querying.
#'
#' @details
#' Sample tables (\code{colTables}) are stored in
#' \code{int_colData(x)$colTables} and provide a mechanism to unnest complex
#' relational data from \code{colData()}. This is particularly useful when
#' sample metadata contains multi-valued properties or hierarchical structures
#' that are difficult to query when embedded as nested DataFrames.
#'
#' Common use cases include:
#' \itemize{
#' \item Patient metadata: multiple disease diagnoses, treatment histories, or
#'   demographic details per sample
#' \item Experimental metadata: replicate-level measurements, batch details, or
#'   quality control metrics
#' \item Cell lineage: developmental trajectories or cell state transitions
#' }
#'
#' Each table is a \linkS4class{DataFrame} with \code{ncol(x)} rows, maintaining
#' alignment with samples. Tables are automatically subset when the parent
#' SingleCellExperiment is subset by columns.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{colTable(x, type, withDimnames=TRUE)}:}{
#' Retrieves a DataFrame containing sample-level relational data. \code{type}
#' is either a string specifying the name of the table in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired table.
#'
#' If \code{withDimnames=TRUE}, row names of the output DataFrame are replaced
#' with the column names of \code{x}.
#' }
#' \item{\code{colTableNames(x)}:}{
#' Returns a character vector containing the names of all tables in \code{x}.
#' This is guaranteed to be of the same length as the number of tables, though
#' the names may not be unique.
#' }
#' \item{\code{colTables(x, withDimnames=TRUE)}:}{
#' Returns a named \linkS4class{SimpleList} of DataFrames containing one or more
#' tables. Each table has the same number of rows as \code{ncol(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of each DataFrame are replaced with
#' the column names of \code{x}.
#' }
#' }
#'
#' @section Single-table setter:
#' \code{colTable(x, type, withDimnames=TRUE) <- value} will add or replace a
#' table in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the table is added or replaced:
#' \itemize{
#' \item If \code{type} is a numeric scalar, it must be within the range of
#'   existing tables, and \code{value} will be assigned to the table at that
#'   index.
#' \item If \code{type} is a string and a table exists with this name,
#'   \code{value} is assigned to that table. Otherwise a new table with this
#'   name is appended to the existing list of tables.
#' }
#'
#' \code{value} is expected to be a DataFrame or data.frame-like object with
#' number of rows equal to \code{ncol(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of \code{value} are set to
#' \code{colnames(x)}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{colTables(x, withDimnames=TRUE) <- value}:}{
#' Replaces all tables in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of DataFrames
#' or data.frame-like objects with number of rows equal to \code{ncol(x)}.
#'
#' If \code{value} is named, those names will be used to name the tables in
#' \code{x}.
#'
#' If \code{value} is a \linkS4class{Annotated} object, any
#' \code{\link{metadata}} will be retained in \code{colTables(x)}.
#' If \code{value} is a \linkS4class{Vector} object, any \code{\link{mcols}}
#' will also be retained.
#'
#' If \code{withDimnames=TRUE}, row names in each entry of \code{value} are set
#' to \code{colnames(x)}.
#' }
#' \item{\code{colTableNames(x) <- value}:}{
#' Replaces all names for tables in \code{x} with a character vector
#' \code{value}. This should be of length equal to the number of tables
#' currently in \code{x}.
#' }
#' }
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param type String or integer scalar specifying the name or index of the
#'   table to get or set.
#' @param withDimnames Logical scalar indicating whether row names should be
#'   extracted from (getters) or set to (setters) the column names of \code{x}.
#' @param ... Additional arguments, currently ignored.
#' @param value For the getter, a DataFrame-like object with number of rows
#'   equal to \code{ncol(x)}, containing the table data. For \code{colTables<-},
#'   a list of such DataFrames. For \code{colTableNames<-}, a character vector
#'   of names.
#'
#' @return
#' For \code{colTable}, a DataFrame containing sample-level relational data.
#'
#' For \code{colTables}, a \linkS4class{SimpleList} of such DataFrames.
#'
#' For \code{colTableNames}, a character vector of table names.
#'
#' For all setters, \code{x} is returned with the modified tables or names.
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' # Create example SCE
#' sce <- SingleCellExperiment(
#'     assays = list(counts = matrix(rpois(1000, 5), 100, 10))
#' )
#' rownames(sce) <- paste0("Gene", 1:100)
#' colnames(sce) <- paste0("Cell", 1:10)
#'
#' # Add a table of patient disease history (multiple diseases per patient)
#' diseases <- DataFrame(
#'     disease_id = paste0("MONDO:", sample(100000:999999, 10)),
#'     disease_label = paste0("Disease_", 1:10),
#'     onset_age = sample(20:80, 10, replace = TRUE)
#' )
#' colTable(sce, "diseases") <- diseases
#'
#' # Add a table of per-sample QC metrics
#' qc_metrics <- DataFrame(
#'     metric_name = rep(c("n_genes", "n_counts", "pct_mito"), length.out = 10),
#'     value = runif(10, 1000, 5000),
#'     passed = sample(c(TRUE, FALSE), 10, replace = TRUE)
#' )
#' colTable(sce, "qc") <- qc_metrics
#'
#' # Retrieve all tables
#' all_tables <- colTables(sce)
#' names(all_tables)
#'
#' # Retrieve specific table by name
#' dis <- colTable(sce, "diseases")
#' dim(dis)  # 10 samples × 3 columns
#'
#' # Retrieve by index
#' qc <- colTable(sce, 2)
#'
#' # Get table names
#' colTableNames(sce)
#'
#' # Subset SCE - tables are automatically subset
#' sce_sub <- sce[, 1:5]
#' dim(colTable(sce_sub, "diseases"))  # 5 samples × 3 columns
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{colData}} for simple sample metadata.
#'
#' \code{\link{rowTables}} for feature-level nested tables.
#'
#' \code{\link[SingleCellExperiment]{int_colData}} for the internal column
#' metadata storage.
#'
#' @aliases
#' colTable
#' colTable,SingleCellExperiment,missing-method
#' colTable,SingleCellExperiment,numeric-method
#' colTable,SingleCellExperiment,character-method
#' colTables
#' colTables,SingleCellExperiment-method
#' colTableNames
#' colTableNames,SingleCellExperiment-method
#' colTable<-
#' colTable<-,SingleCellExperiment,missing-method
#' colTable<-,SingleCellExperiment,numeric-method
#' colTable<-,SingleCellExperiment,character-method
#' colTables<-
#' colTables<-,SingleCellExperiment-method
#' colTableNames<-
#' colTableNames<-,SingleCellExperiment,character-method
#'
#' @include SingleCellExperiments-internals.R
#'
#' @name SingleCellExperiment-colTables
NULL

.ctabs_key <- "colTables"

#' @importClassesFrom S4Vectors DataFrame
.any2dataframe <- function(value) {
    if (!is(value, "DataFrame")) {
        value <- as(value, "DataFrame")
    }
    value
}

#' @importFrom S4Vectors make_zero_col_DFrame
.initialize_colTables <- function(x) {
    if (is.null(x@int_colData[[.ctabs_key]])) {
        tables <- make_zero_col_DFrame(ncol(x))
        rownames(tables) <- colnames(x)
        x@int_colData[[.ctabs_key]] <- tables
    }
    x
}

#' @export
setGeneric("colTables", function(x, withDimnames=TRUE) standardGeneric("colTables"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("colTables", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    x <- .initialize_colTables(x)
    value <- .get_internal_all(x, getfun = int_colData, key = .ctabs_key)
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- colnames(x)
        }
    }
    value
})

#' @export
setGeneric("colTables<-", function(x, withDimnames=TRUE, ..., value) standardGeneric("colTables<-"))

#' @export
#' @importClassesFrom S4Vectors List
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("colTables", "SingleCellExperiment", function(x, withDimnames=TRUE, ..., value) {
    if (!(is.null(value) || is.list(value) || is(value, "List"))) {
        stop("'invalid 'value' in 'colTables(<SingleCellExperiment>) <- value'")
    }

    for (i in seq_along(value)) {
        if (!identical(nrow(value[[i]]), ncol(x))) {
            stop("invalid 'value' in 'colTables(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "all elements in 'value' should have number of rows equal to 'nrow(x)'")
        }
        if (withDimnames && !is.null(rownames(value[[i]]))) {
            if (!identical(rownames(value[[i]]), colnames(x))) {
                warning("non-NULL 'rownames(value[[", i, "]])' should be the same as 'colnames(x)' for 'colTables<-'.")
            }
        }
    }

    x <- .initialize_colTables(x)
    .set_internal_all(x, value,
                      getfun = int_colData,
                      setfun = `int_colData<-`,
                      key = .ctabs_key,
                      convertfun = .any2dataframe,
                      xdimfun = ncol,
                      vdimfun = nrow,
                      funstr = "colTables",
                      xdimstr = "ncol",
                      vdimstr = "rows")
})

#' @export
setGeneric("colTableNames", function(x) standardGeneric("colTableNames"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("colTableNames", "SingleCellExperiment", function(x) {
    x <- .initialize_colTables(x)
    .get_internal_names(x, getfun = int_colData, key = .ctabs_key)
})

#' @export
setGeneric("colTableNames<-", function(x, value) standardGeneric("colTableNames<-"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("colTableNames", c("SingleCellExperiment", "character"), function(x, value) {
    x <- .initialize_colTables(x)
    .set_internal_names(x, value,
                        getfun = int_colData,
                        setfun = `int_colData<-`,
                        key = .ctabs_key)
})

#' @export
setGeneric("colTable", function(x, type, withDimnames=TRUE) standardGeneric("colTable"))

#' @export
setMethod("colTable", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_colTables(x)
    .get_internal_missing(x,
                          basefun = colTable,
                          namefun = colTableNames,
                          funstr = "colTable",
                          withDimnames = withDimnames)
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("colTable", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_colTables(x)
    out <- .get_internal_integer(x, type,
                                 getfun = int_colData,
                                 key = .ctabs_key,
                                 funstr = "colTable",
                                 substr = "type")
    if (withDimnames) {
        rownames(out) <- colnames(x)
    }
    out
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("colTable", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_colTables(x)
    out <- .get_internal_character(x, type,
                                   getfun = int_colData,
                                   key = .ctabs_key,
                                   funstr = "colTable",
                                   substr = "type",
                                   namestr = "colTableNames")
    if (withDimnames) {
        rownames(out) <- colnames(x)
    }
    out
})

#' @export
setGeneric("colTable<-", function(x, type, withDimnames=TRUE, ..., value) standardGeneric("colTable<-"))

#' @export
setReplaceMethod("colTable", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE, ..., value) {
    x <- .initialize_colTables(x)
    .set_internal_missing(x, value,
                          withDimnames = withDimnames,
                          basefun = `colTable<-`,
                          namefun = colTableNames)
})


#' @export
#' @importFrom S4Vectors isSingleNumber
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("colTable", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleNumber(type)) {
        stop("'type' must be a scalar in 'colTable(<SingleCellExperiment>, type=\"numeric\")' <- value'")
    }

    if (type > length(colTableNames(x))) {
        stop("'type' out of bounds in 'colTable(<SingleCellExperiment>, type=\"numeric\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), ncol(x))) {
            stop("invalid 'value' in 'colTable(<SingleCellExperiment>, type=\"numeric\") <- value':\n  ",
                 "'value' should have number of rows equal to 'ncol(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), colnames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'colnames(x)' for 'colTable<-'.")
            }
        }
    }

    x <- .initialize_colTables(x)
    .set_internal_numeric(x, type, value,
                          getfun = int_colData,
                          setfun = `int_colData<-`,
                          key = .ctabs_key,
                          convertfun = .any2dataframe,
                          xdimfun = ncol,
                          vdimfun = nrow,
                          funstr = "colTable",
                          xdimstr = "ncol",
                          vdimstr = "rows",
                          substr = "type")
})

#' @export
#' @importFrom S4Vectors isSingleString
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("colTable", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleString(type)) {
        stop("'type' must be a string in 'colTable(<SingleCellExperiment>, type=\"character\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), ncol(x))) {
            stop("invalid 'value' in 'colTable(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "'value' should have number of rows equal to 'ncol(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), colnames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'colnames(x)' for 'colTable<-'.")
            }
        }
    }

    x <- .initialize_colTables(x)
    .set_internal_character(x, type, value,
                            getfun = int_colData,
                            setfun = `int_colData<-`,
                            key = .ctabs_key,
                            convertfun = .any2dataframe,
                            xdimfun = ncol,
                            vdimfun = nrow,
                            funstr = "colTable",
                            xdimstr = "ncol",
                            vdimstr = "rows",
                            substr = "type")
})
