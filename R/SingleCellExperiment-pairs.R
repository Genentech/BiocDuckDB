#' colPairs/rowPairs setters for DuckDBSelfHits
#'
#' @description
#' Specialized setter methods for \code{colPairs} and \code{rowPairs} that
#' automatically wrap \linkS4class{DuckDBSelfHits} objects in
#' \linkS4class{DuckDBDualSubset} containers. This enables lazy node-based
#' subsetting when a SingleCellExperiment is subset.
#'
#' @details
#' These methods extend the default \code{colPairs<-} and \code{rowPairs<-}
#' setters from SingleCellExperiment to recognize \code{DuckDBSelfHits} objects
#' and wrap them in \code{DuckDBDualSubset} containers. This ensures that when
#' a SingleCellExperiment is subset (e.g., \code{sce[, 1:100]}), the wrapped
#' \code{DuckDBSelfHits} objects are lazily filtered via SQL rather than
#' materialized in memory.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param type String or integer scalar specifying the name or index of the
#'   pairing to replace.
#' @param ... Additional arguments passed to the default setter.
#' @param value A \linkS4class{DuckDBSelfHits} object to be stored as a
#'   column/row pairing. This will be automatically wrapped in a
#'   \linkS4class{DuckDBDualSubset}.
#'
#' @return
#' For the setter methods, \code{x} is returned with the modified pairings.
#'
#' @section Usage:
#' \preformatted{
#' # Create SCE with lazy KNN graph
#' library(SingleCellExperiment)
#' library(BiocDuckDB)
#'
#' sce <- SingleCellExperiment(assays = list(counts = matrix(1:100, 10, 10)))
#'
#' # Create DuckDBSelfHits from parquet
#' knn_hits <- DuckDBSelfHits("knn_edges.parquet",
#'                            from = "from",
#'                            to = "to",
#'                            nnode = ncol(sce))
#'
#' # Assign to colPairs (auto-wraps in DuckDBDualSubset)
#' colPairs(sce, "knn") <- knn_hits
#'
#' # Subsetting SCE automatically subsets the graph lazily
#' sce_sub <- sce[, 1:5]
#' pairs_sub <- colPairs(sce_sub, "knn")  # Lazy filtered to nodes 1:5
#' }
#'
#' @author Patrick Aboyoun
#'
#' @aliases
#' colPair<-,SingleCellExperiment,character,DuckDBSelfHits-method
#' colPair<-,SingleCellExperiment,missing,DuckDBSelfHits-method
#' colPair<-,SingleCellExperiment,numeric,DuckDBSelfHits-method
#' colPairs<-,SingleCellExperiment,DuckDBSelfHits-method
#' rowPair<-,SingleCellExperiment,character,DuckDBSelfHits-method
#' rowPair<-,SingleCellExperiment,missing,DuckDBSelfHits-method
#' rowPair<-,SingleCellExperiment,numeric,DuckDBSelfHits-method
#' rowPairs<-,SingleCellExperiment,DuckDBSelfHits-method
#'
#' @seealso
#' \linkS4class{DuckDBSelfHits} for the lazy edge storage class.
#'
#' \linkS4class{DuckDBDualSubset} for the wrapper that enables node-based subsetting.
#'
#' \code{\link[SingleCellExperiment]{colPairs}} and
#' \code{\link[SingleCellExperiment]{rowPairs}} for the base SingleCellExperiment methods.
#'
#' @name SingleCellExperiment-pairs
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### colPairs setters
###

#' @export
#' @importFrom SingleCellExperiment colPair<- colPairs<-
setReplaceMethod("colPairs", c("SingleCellExperiment", "DuckDBSelfHits"),
function(x, value) {
    colPair(x) <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment colPair<- colPairNames
setReplaceMethod("colPair", c("SingleCellExperiment", "missing", "DuckDBSelfHits"),
function(x, type, ..., value) {
    type <- if (length(colPairNames(x))) 1L else "unnamed1"
    colPair(x, type) <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment colPair<-
#' @importFrom S4Vectors make_zero_col_DFrame
setReplaceMethod("colPair", c("SingleCellExperiment", "numeric", "DuckDBSelfHits"),
function(x, type, ..., value) {
    value <- DuckDBDualSubset(value)
    if (is.null(x@int_colData@listData[["colPairs"]])) {
        df <- make_zero_col_DFrame(length(value))
        x@int_colData@listData[["colPairs"]] <- df
    }
    x@int_colData@listData[["colPairs"]][[type]] <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment colPair<-
#' @importFrom S4Vectors make_zero_col_DFrame
setReplaceMethod("colPair", c("SingleCellExperiment", "character", "DuckDBSelfHits"),
function(x, type, ..., value) {
    value <- DuckDBDualSubset(value)
    if (is.null(x@int_colData@listData[["colPairs"]])) {
        df <- make_zero_col_DFrame(length(value))
        x@int_colData@listData[["colPairs"]] <- df
    }
    x@int_colData@listData[["colPairs"]][[type]] <- value
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### rowPairs setters
###

#' @export
#' @importFrom SingleCellExperiment rowPair<- rowPairs<-
setReplaceMethod("rowPairs", c("SingleCellExperiment", "DuckDBSelfHits"),
function(x, value) {
    rowPair(x) <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment rowPair<- rowPairNames
setReplaceMethod("rowPair", c("SingleCellExperiment", "missing", "DuckDBSelfHits"),
function(x, type, ..., value) {
    type <- if (length(rowPairNames(x))) 1L else "unnamed1"
    rowPair(x, type) <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment rowPair<-
#' @importFrom S4Vectors make_zero_col_DFrame
setReplaceMethod("rowPair", c("SingleCellExperiment", "numeric", "DuckDBSelfHits"),
function(x, type, ..., value) {
    value <- DuckDBDualSubset(value)
    if (is.null(x@int_elementMetadata@listData[["rowPairs"]])) {
        df <- make_zero_col_DFrame(length(value))
        x@int_elementMetadata@listData[["rowPairs"]] <- df
    }
    x@int_elementMetadata@listData[["rowPairs"]][[type]] <- value
    x
})

#' @export
#' @importFrom SingleCellExperiment rowPair<-
#' @importFrom S4Vectors make_zero_col_DFrame
setReplaceMethod("rowPair", c("SingleCellExperiment", "character", "DuckDBSelfHits"),
function(x, type, ..., value) {
    value <- DuckDBDualSubset(value)
    if (is.null(x@int_elementMetadata@listData[["rowPairs"]])) {
        df <- make_zero_col_DFrame(length(value))
        x@int_elementMetadata@listData[["rowPairs"]] <- df
    }
    x@int_elementMetadata@listData[["rowPairs"]][[type]] <- value
    x
})
