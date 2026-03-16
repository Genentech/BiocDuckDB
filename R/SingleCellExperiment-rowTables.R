#' Feature table methods
#'
#' @description
#' Methods to get or set nested feature-level tables in a
#' \linkS4class{SingleCellExperiment} object. These tables enable storage of
#' multi-column relational data associated with features (genes), which would
#' otherwise require nested DataFrames in \code{rowData()}. By extracting nested
#' tables into a separate slot, they can be serialized to independent Parquet
#' tables for efficient querying.
#'
#' @details
#' Feature tables (\code{rowTables}) are stored in
#' \code{int_elementMetadata(x)$rowTables} and provide a mechanism to unnest
#' complex relational data from \code{rowData()}. This is particularly useful
#' when feature metadata contains multi-valued properties or hierarchical
#' structures that are difficult to query when embedded as nested DataFrames.
#'
#' Common use cases include:
#' \itemize{
#' \item Gene annotation tables: transcript isoforms, protein domains, or GO terms
#' \item Experimental design tables: per-gene treatment conditions or batch effects
#' \item Statistical results: detailed per-gene statistics across multiple contrasts
#' }
#'
#' Each table is a \linkS4class{DataFrame} with \code{nrow(x)} rows, maintaining
#' alignment with features. Tables are automatically subset when the parent
#' SingleCellExperiment is subset by rows.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{rowTable(x, type, withDimnames=TRUE)}:}{
#' Retrieves a DataFrame containing feature-level relational data. \code{type}
#' is either a string specifying the name of the table in \code{x} to retrieve,
#' or a numeric scalar specifying the index of the desired table.
#'
#' If \code{withDimnames=TRUE}, row names of the output DataFrame are replaced
#' with the row names of \code{x}.
#' }
#' \item{\code{rowTableNames(x)}:}{
#' Returns a character vector containing the names of all tables in \code{x}.
#' This is guaranteed to be of the same length as the number of tables, though
#' the names may not be unique.
#' }
#' \item{\code{rowTables(x, withDimnames=TRUE)}:}{
#' Returns a named \linkS4class{SimpleList} of DataFrames containing one or more
#' tables. Each table has the same number of rows as \code{nrow(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of each DataFrame are replaced with
#' the row names of \code{x}.
#' }
#' }
#'
#' @section Single-table setter:
#' \code{rowTable(x, type, withDimnames=TRUE) <- value} will add or replace a
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
#' number of rows equal to \code{nrow(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of \code{value} are set to
#' \code{rownames(x)}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{rowTables(x, withDimnames=TRUE) <- value}:}{
#' Replaces all tables in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of DataFrames
#' or data.frame-like objects with number of rows equal to \code{nrow(x)}.
#'
#' If \code{value} is named, those names will be used to name the tables in
#' \code{x}.
#'
#' If \code{value} is a \linkS4class{Annotated} object, any
#' \code{\link{metadata}} will be retained in \code{rowTables(x)}.
#' If \code{value} is a \linkS4class{Vector} object, any \code{\link{mcols}}
#' will also be retained.
#'
#' If \code{withDimnames=TRUE}, row names in each entry of \code{value} are set
#' to \code{rownames(x)}.
#' }
#' \item{\code{rowTableNames(x) <- value}:}{
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
#'   extracted from (getters) or set to (setters) the row names of \code{x}.
#' @param ... Additional arguments, currently ignored.
#' @param value For the getter, a DataFrame-like object with number of rows
#'   equal to \code{nrow(x)}, containing the table data. For \code{rowTables<-},
#'   a list of such DataFrames. For \code{rowTableNames<-}, a character vector
#'   of names.
#'
#' @return
#' For \code{rowTable}, a DataFrame containing feature-level relational data.
#'
#' For \code{rowTables}, a \linkS4class{SimpleList} of such DataFrames.
#'
#' For \code{rowTableNames}, a character vector of table names.
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
#' # Add a table of gene isoforms (multiple isoforms per gene)
#' isoforms <- DataFrame(
#'     isoform_id = paste0("ISO", 1:100),
#'     length = sample(500:5000, 100),
#'     is_canonical = sample(c(TRUE, FALSE), 100, replace = TRUE)
#' )
#' rowTable(sce, "isoforms") <- isoforms
#'
#' # Add a table of differential expression results across conditions
#' de_results <- DataFrame(
#'     condition = rep(c("treated", "control"), each = 50),
#'     log2fc = rnorm(100),
#'     pvalue = runif(100)
#' )
#' rowTable(sce, "de_stats") <- de_results
#'
#' # Retrieve all tables
#' all_tables <- rowTables(sce)
#' names(all_tables)
#'
#' # Retrieve specific table by name
#' iso <- rowTable(sce, "isoforms")
#' dim(iso)  # 100 genes × 3 columns
#'
#' # Retrieve by index
#' de <- rowTable(sce, 2)
#'
#' # Get table names
#' rowTableNames(sce)
#'
#' # Subset SCE - tables are automatically subset
#' sce_sub <- sce[1:50, ]
#' dim(rowTable(sce_sub, "isoforms"))  # 50 genes × 3 columns
#'
#' @seealso
#' \code{\link[SummarizedExperiment]{rowData}} for simple feature metadata.
#'
#' \code{\link{colTables}} for sample-level nested tables.
#'
#' \code{\link[SingleCellExperiment]{int_elementMetadata}} for the internal row
#' metadata storage.
#'
#' @aliases
#' rowTable
#' rowTable,SingleCellExperiment,missing-method
#' rowTable,SingleCellExperiment,numeric-method
#' rowTable,SingleCellExperiment,character-method
#' rowTables
#' rowTables,SingleCellExperiment-method
#' rowTableNames
#' rowTableNames,SingleCellExperiment-method
#' rowTable<-
#' rowTable<-,SingleCellExperiment,missing-method
#' rowTable<-,SingleCellExperiment,numeric-method
#' rowTable<-,SingleCellExperiment,character-method
#' rowTables<-
#' rowTables<-,SingleCellExperiment-method
#' rowTableNames<-
#' rowTableNames<-,SingleCellExperiment,character-method
#'
#' @include SingleCellExperiments-internals.R
#'
#' @name SingleCellExperiment-rowTables
NULL

.rtabs_key <- "rowTables"

#' @importClassesFrom S4Vectors DataFrame
.any2dataframe <- function(value) {
    if (!is(value, "DataFrame")) {
        value <- as(value, "DataFrame")
    }
    value
}

#' @importFrom S4Vectors make_zero_col_DFrame
.initialize_rowTables <- function(x) {
    if (is.null(x@int_elementMetadata[[.rtabs_key]])) {
        tables <- make_zero_col_DFrame(nrow(x))
        rownames(tables) <- rownames(x)
        x@int_elementMetadata[[.rtabs_key]] <- tables
    }
    x
}

#' @export
setGeneric("rowTables", function(x, withDimnames=TRUE) standardGeneric("rowTables"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowTables", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    x <- .initialize_rowTables(x)
    value <- .get_internal_all(x, getfun = int_elementMetadata, key = .rtabs_key)
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- rownames(x)
        }
    }
    value
})

#' @export
setGeneric("rowTables<-", function(x, withDimnames=TRUE, ..., value) standardGeneric("rowTables<-"))

#' @export
#' @importClassesFrom S4Vectors List
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowTables", "SingleCellExperiment", function(x, withDimnames=TRUE, ..., value) {
    if (!(is.null(value) || is.list(value) || is(value, "List"))) {
        stop("'invalid 'value' in 'rowTables(<SingleCellExperiment>) <- value'")
    }

    for (i in seq_along(value)) {
        if (!identical(nrow(value[[i]]), nrow(x))) {
            stop("invalid 'value' in 'rowTables(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "all elements in 'value' should have number of rows equal to 'nrow(x)'")
        }
        if (withDimnames && !is.null(rownames(value[[i]]))) {
            if (!identical(rownames(value[[i]]), rownames(x))) {
                warning("non-NULL 'rownames(value[[", i, "]])' should be the same as 'rownames(x)' for 'rowTables<-'.")
            }
        }
    }

    x <- .initialize_rowTables(x)
    .set_internal_all(x, value,
                      getfun = int_elementMetadata,
                      setfun = `int_elementMetadata<-`,
                      key = .rtabs_key,
                      convertfun = .any2dataframe,
                      xdimfun = nrow,
                      vdimfun = nrow,
                      funstr = "rowTables",
                      xdimstr = "nrow",
                      vdimstr = "rows")
})

#' @export
setGeneric("rowTableNames", function(x) standardGeneric("rowTableNames"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowTableNames", "SingleCellExperiment", function(x) {
    x <- .initialize_rowTables(x)
    .get_internal_names(x, getfun = int_elementMetadata, key = .rtabs_key)
})

#' @export
setGeneric("rowTableNames<-", function(x, value) standardGeneric("rowTableNames<-"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowTableNames", c("SingleCellExperiment", "character"), function(x, value) {
    x <- .initialize_rowTables(x)
    .set_internal_names(x, value,
                        getfun = int_elementMetadata,
                        setfun = `int_elementMetadata<-`,
                        key = .rtabs_key)
})

#' @export
setGeneric("rowTable", function(x, type, withDimnames=TRUE) standardGeneric("rowTable"))

#' @export
setMethod("rowTable", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_rowTables(x)
    .get_internal_missing(x,
                          basefun = rowTable,
                          namefun = rowTableNames,
                          funstr = "rowTable",
                          withDimnames = withDimnames)
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowTable", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_rowTables(x)
    out <- .get_internal_integer(x, type,
                                 getfun = int_elementMetadata,
                                 key = .rtabs_key,
                                 funstr = "rowTable",
                                 substr = "type")
    if (withDimnames) {
        rownames(out) <- rownames(x)
    }
    out
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowTable", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_rowTables(x)
    out <- .get_internal_character(x, type,
                                   getfun = int_elementMetadata,
                                   key = .rtabs_key,
                                   funstr = "rowTable",
                                   substr = "type",
                                   namestr = "rowTableNames")
    if (withDimnames) {
        rownames(out) <- rownames(x)
    }
    out
})

#' @export
setGeneric("rowTable<-", function(x, type, withDimnames=TRUE, ..., value) standardGeneric("rowTable<-"))

#' @export
setReplaceMethod("rowTable", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE, ..., value) {
    x <- .initialize_rowTables(x)
    .set_internal_missing(x, value,
                          withDimnames = withDimnames,
                          basefun = `rowTable<-`,
                          namefun = rowTableNames)
})


#' @export
#' @importFrom S4Vectors isSingleNumber
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowTable", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleNumber(type)) {
        stop("'type' must be a scalar in 'rowTable(<SingleCellExperiment>, type=\"numeric\")' <- value'")
    }

    if (type > length(rowTableNames(x))) {
        stop("'type' out of bounds in 'rowTable(<SingleCellExperiment>, type=\"numeric\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), nrow(x))) {
            stop("invalid 'value' in 'rowTable(<SingleCellExperiment>, type=\"numeric\") <- value':\n  ",
                 "'value' should have number of rows equal to 'nrow(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), rownames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowTable<-'.")
            }
        }
    }

    x <- .initialize_rowTables(x)
    .set_internal_numeric(x, type, value,
                          getfun = int_elementMetadata,
                          setfun = `int_elementMetadata<-`,
                          key = .rtabs_key,
                          convertfun = .any2dataframe,
                          xdimfun = nrow,
                          vdimfun = nrow,
                          funstr = "rowTable",
                          xdimstr = "nrow",
                          vdimstr = "rows",
                          substr = "type")
})

#' @export
#' @importFrom S4Vectors isSingleString
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowTable", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleString(type)) {
        stop("'type' must be a string in 'rowTable(<SingleCellExperiment>, type=\"character\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), nrow(x))) {
            stop("invalid 'value' in 'rowTable(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "'value' should have number of rows equal to 'nrow(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), rownames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowTable<-'.")
            }
        }
    }

    x <- .initialize_rowTables(x)
    .set_internal_character(x, type, value,
                            getfun = int_elementMetadata,
                            setfun = `int_elementMetadata<-`,
                            key = .rtabs_key,
                            convertfun = .any2dataframe,
                            xdimfun = nrow,
                            vdimfun = nrow,
                            funstr = "rowTable",
                            xdimstr = "nrow",
                            vdimstr = "rows",
                            substr = "type")
})
