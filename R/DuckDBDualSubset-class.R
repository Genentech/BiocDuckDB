#' DuckDBDualSubset objects
#'
#' @description
#' The DuckDBDualSubset class conforms to the DualSubset class from
#' SingleCellExperiment to provide lazy node-based subsetting for
#' \linkS4class{DuckDBSelfHits} objects.
#'
#' @details
#' DuckDBDualSubset is a lightweight wrapper around \linkS4class{DuckDBSelfHits} 
#' that provides node-based subsetting semantics required by SingleCellExperiment's
#' \code{colPairs}/\code{rowPairs} framework. Unlike the \code{DualSubset} class
#' which materializes data immediately upon subsetting, \code{DuckDBDualSubset}
#' delegates to the \code{nodes} slot in the wrapped \code{DuckDBSelfHits} object, 
#' enabling lazy SQL-based edge filtering.
#'
#' When a \code{DuckDBDualSubset} is subset via \code{x[i]}, it delegates to
#' \code{extractNODES(x@@hits, i)}, which updates the \code{nodes} slot in the
#' \code{DuckDBSelfHits} object. Edge filtering is deferred until materialization
#' via SQL WHERE clauses (only edges where BOTH from and to are in the node subset
#' are retained).
#'
#' @section Constructor:
#' \describe{
#'   \item{\code{DuckDBDualSubset(hits)}:}{
#'     Creates a DuckDBDualSubset object wrapping a DuckDBSelfHits object.
#'     \describe{
#'       \item{\code{hits}}{
#'         A \linkS4class{DuckDBSelfHits} object containing the edge data.
#'       }
#'     }
#'   }
#' }
#'
#' @section Accessors:
#' \describe{
#'   \item{\code{length(x)}:}{
#'     Returns the number of nodes (\code{nnode(x@@hits)}).
#'   }
#' }
#'
#' @section Subsetting:
#' \describe{
#'   \item{\code{x[i]}:}{
#'     Performs node-based subsetting by delegating to
#'     \code{extractNODES(x@@hits, i)}. This updates the \code{nodes} slot
#'     in the wrapped \code{DuckDBSelfHits} object and applies lazy SQL filtering
#'     to retain only edges where BOTH endpoints are in the subset.
#'   }
#' }
#'
#' @section Usage with SingleCellExperiment:
#' DuckDBDualSubset objects are typically created automatically when assigning
#' DuckDBSelfHits objects to \code{colPairs}/\code{rowPairs}:
#'
#' @author Patrick Aboyoun
#'
#' @aliases [,DuckDBDualSubset,ANY,ANY,ANY-method
#' @aliases [<-,DuckDBDualSubset,ANY,ANY,ANY-method
#' @aliases c,DuckDBDualSubset-method
#' @aliases length,DuckDBDualSubset-method
#'
#' @seealso
#' \linkS4class{DuckDBSelfHits} for the underlying edge storage class.
#'
#' \code{\link[DuckDBDataFrame]{extractNODES}} for the node-based subsetting method.
#'
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(counts = matrix(1:50, 10, 5)))
#' hits <- S4Vectors::SelfHits(from = 1:3, to = 2:4, nnode = 5L)
#' tmp <- tempfile()
#' writeParquet(hits, tmp)
#' ddb_hits <- DuckDBSelfHits(tmp, from = "from", to = "to", nnode = 5L)
#' colPair(sce, "knn") <- ddb_hits
#' ds <- colPair(sce, "knn")
#' length(ds)
#' ds[1:3]
#' unlink(tmp, recursive = TRUE)
#'
#' @name DuckDBDualSubset-class
NULL

replaceSlots <- BiocGenerics:::replaceSlots

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Class definition
###

#' @importClassesFrom DuckDBDataFrame DuckDBSelfHits
setClass("DuckDBDualSubset", slots = c(hits = "DuckDBSelfHits"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @importFrom S4Vectors new2
DuckDBDualSubset <- function(hits) {
    if (!is(hits, "DuckDBSelfHits")) {
        stop("'hits' must be a DuckDBSelfHits object")
    }
    new2("DuckDBDualSubset", hits = hits, check = FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

#' @importFrom S4Vectors nnode
setMethod("length", "DuckDBDualSubset", function(x) nnode(x@hits))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

#' @importFrom DuckDBDataFrame extractNODES
setMethod("[", "DuckDBDualSubset", function(x, i, j, ..., drop = FALSE) {
    if (!missing(j) || length(list(...)) > 0L) {
        stop("DuckDBDualSubset subsetting only supports 'x[i]'")
    }
    if (missing(i)) {
        return(x)
    }
    replaceSlots(x, hits = extractNODES(x@hits, i), check = FALSE)
})

setReplaceMethod("[", "DuckDBDualSubset", function(x, i, j, ..., value) {
    stop("element replacement is not supported for DuckDBDualSubset objects")
})

setMethod("c", "DuckDBDualSubset", function(x, ...) {
    stop("concatenation of DuckDBDualSubset objects is not supported")
})
