#' Feature loading methods
#'
#' @description
#' Methods to get or set feature-level loading matrices in a
#' \linkS4class{SingleCellExperiment} object. These matrices represent
#' per-feature contributions to latent factors, embeddings, or statistical
#' models, where feature identity (e.g., gene names) carries semantic meaning.
#' Each row of a loading matrix corresponds to a row (feature) of the
#' SingleCellExperiment object.
#'
#' @details
#' Feature loadings (\code{rowLoadings}) are stored in
#' \code{int_elementMetadata(x)$rowLoadings} and represent the AnnData
#' \code{.varm} equivalent in R. Unlike sample embeddings (\code{reducedDims})
#' where dimension labels are arbitrary coordinates, loading matrices have
#' \strong{meaningful row names} that identify which features contribute to each
#' latent dimension.
#'
#' Common use cases include:
#' \itemize{
#' \item PCA loadings: which genes load on each principal component
#' \item Factor analysis: gene weights for latent factors
#' \item Statistical summaries: per-gene effect sizes or log fold changes across
#'   conditions
#' }
#'
#' Loading matrices are automatically subset when the parent
#' SingleCellExperiment is subset by rows, maintaining alignment between
#' features and their loadings.
#'
#' @section Getters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{rowLoading(x, type, withDimnames=TRUE)}:}{
#' Retrieves a matrix containing feature loadings (rows) for each latent
#' dimension (columns). \code{type} is either a string specifying the name of
#' the loading result in \code{x} to retrieve, or a numeric scalar specifying
#' the index of the desired result.
#'
#' If \code{withDimnames=TRUE}, row names of the output matrix are replaced with
#' the row names of \code{x}.
#' }
#' \item{\code{rowLoadingNames(x)}:}{
#' Returns a character vector containing the names of all loading results in
#' \code{x}. This is guaranteed to be of the same length as the number of
#' results, though the names may not be unique.
#' }
#' \item{\code{rowLoadings(x, withDimnames=TRUE)}:}{
#' Returns a named \linkS4class{SimpleList} of matrices containing one or more
#' loading results. Each result is a matrix with the same number of rows as
#' \code{nrow(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of each matrix are replaced with the
#' row names of \code{x}.
#' }
#' }
#'
#' @section Single-result setter:
#' \code{rowLoading(x, type, withDimnames=TRUE) <- value} will add or replace a
#' loading result in a \linkS4class{SingleCellExperiment} object \code{x}.
#' The value of \code{type} determines how the result is added or replaced:
#' \itemize{
#' \item If \code{type} is a numeric scalar, it must be within the range of
#'   existing results, and \code{value} will be assigned to the result at that
#'   index.
#' \item If \code{type} is a string and a result exists with this name,
#'   \code{value} is assigned to that result. Otherwise a new result with this
#'   name is appended to the existing list of results.
#' }
#'
#' \code{value} is expected to be a matrix or matrix-like object with number of
#' rows equal to \code{nrow(x)}.
#'
#' If \code{withDimnames=TRUE}, row names of \code{value} are set to
#' \code{rownames(x)}.
#'
#' @section Other setters:
#' In the following examples, \code{x} is a \linkS4class{SingleCellExperiment}
#' object.
#' \describe{
#' \item{\code{rowLoadings(x, withDimnames=TRUE) <- value}:}{
#' Replaces all loading results in \code{x} with those in \code{value}.
#' The latter should be a list-like object containing any number of matrices or
#' matrix-like objects with number of rows equal to \code{nrow(x)}.
#'
#' If \code{value} is named, those names will be used to name the loading
#' results in \code{x}.
#'
#' If \code{value} is a \linkS4class{Annotated} object, any
#' \code{\link{metadata}} will be retained in \code{rowLoadings(x)}.
#' If \code{value} is a \linkS4class{Vector} object, any \code{\link{mcols}}
#' will also be retained.
#'
#' If \code{withDimnames=TRUE}, row names in each entry of \code{value} are set
#' to \code{rownames(x)}.
#' }
#' \item{\code{rowLoadingNames(x) <- value}:}{
#' Replaces all names for loading results in \code{x} with a character vector
#' \code{value}. This should be of length equal to the number of results
#' currently in \code{x}.
#' }
#' }
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' @param type String or integer scalar specifying the name or index of the
#'   loading result to get or set.
#' @param withDimnames Logical scalar indicating whether row names should be
#'   extracted from (getters) or set to (setters) the row names of \code{x}.
#' @param ... Additional arguments, currently ignored.
#' @param value For the getter, a matrix-like object with number of rows equal
#'   to \code{nrow(x)}, containing the loading values. For \code{rowLoadings<-},
#'   a list of such matrices. For \code{rowLoadingNames<-}, a character vector
#'   of names.
#'
#' @return
#' For \code{rowLoading}, a matrix containing loading values for features (rows)
#' and dimensions (columns).
#'
#' For \code{rowLoadings}, a \linkS4class{SimpleList} of such matrices.
#'
#' For \code{rowLoadingNames}, a character vector of loading names.
#'
#' For all setters, \code{x} is returned with the modified loading results or
#' names.
#'
#' @author Patrick Aboyoun
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' # Create example SCE with perturbation data
#' sce <- SingleCellExperiment(
#'     assays = list(beta_Z = matrix(rnorm(1000), 100, 10))
#' )
#' rownames(sce) <- paste0("Gene", 1:100)
#' colnames(sce) <- paste0("Pert", 1:10)
#'
#' # Add PCA loadings (genes × components)
#' pca_loadings <- matrix(rnorm(100 * 5), 100, 5)
#' rownames(pca_loadings) <- rownames(sce)
#' colnames(pca_loadings) <- paste0("PC", 1:5)
#' rowLoading(sce, "PCA") <- pca_loadings
#'
#' # Add factor analysis loadings (genes × factors)
#' factor_loadings <- matrix(rnorm(100 * 20), 100, 20)
#' rownames(factor_loadings) <- rownames(sce)
#' colnames(factor_loadings) <- paste0("Factor", 1:20)
#' rowLoading(sce, "FactorAnalysis") <- factor_loadings
#'
#' # Retrieve all loadings
#' all_loadings <- rowLoadings(sce)
#' names(all_loadings)
#'
#' # Retrieve specific loading by name
#' pca <- rowLoading(sce, "PCA")
#' dim(pca)  # 100 genes × 5 components
#'
#' # Retrieve by index
#' factors <- rowLoading(sce, 2)
#'
#' # Get loading names
#' rowLoadingNames(sce)
#'
#' # Subset SCE - loadings are automatically subset
#' sce_sub <- sce[1:50, ]
#' dim(rowLoading(sce_sub, "PCA"))  # 50 genes × 5 components
#'
#' @seealso
#' \code{\link[SingleCellExperiment]{reducedDims}} for sample-level embeddings
#' (observation matrices).
#'
#' \code{\link[SingleCellExperiment]{int_elementMetadata}} for the internal row
#' metadata storage.
#'
#' @aliases rowLoading
#' @aliases rowLoading,SingleCellExperiment,missing-method
#' @aliases rowLoading,SingleCellExperiment,numeric-method
#' @aliases rowLoading,SingleCellExperiment,character-method
#' @aliases rowLoadings
#' @aliases rowLoadings,SingleCellExperiment-method
#' @aliases rowLoadingNames
#' @aliases rowLoadingNames,SingleCellExperiment-method
#' @aliases rowLoading<-
#' @aliases rowLoading<-,SingleCellExperiment,missing-method
#' @aliases rowLoading<-,SingleCellExperiment,numeric-method
#' @aliases rowLoading<-,SingleCellExperiment,character-method
#' @aliases rowLoadings<-
#' @aliases rowLoadings<-,SingleCellExperiment-method
#' @aliases rowLoadingNames<-
#' @aliases rowLoadingNames<-,SingleCellExperiment,character-method
#'
#' @include SingleCellExperiments-internals.R
#'
#' @name SingleCellExperiment-loadings
NULL

.load_key <- "rowLoadings"

#' @importFrom S4Vectors make_zero_col_DFrame
.initialize_loadings <- function(x) {
    if (is.null(x@int_elementMetadata[[.load_key]])) {
        loadings <- make_zero_col_DFrame(nrow(x))
        rownames(loadings) <- rownames(x)
        x@int_elementMetadata[[.load_key]] <- loadings
    }
    x
}

#' @export
setGeneric("rowLoadings", function(x, withDimnames=TRUE) standardGeneric("rowLoadings"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowLoadings", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    x <- .initialize_loadings(x)
    value <- .get_internal_all(x, getfun = int_elementMetadata, key = .load_key)
    if (withDimnames) {
        for (i in seq_along(value)) {
            rownames(value[[i]]) <- rownames(x)
        }
    }
    value
})

#' @export
setGeneric("rowLoadings<-", function(x, withDimnames=TRUE, ..., value) standardGeneric("rowLoadings<-"))

#' @export
#' @importClassesFrom S4Vectors List
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowLoadings", "SingleCellExperiment", function(x, withDimnames=TRUE, ..., value) {
    if (!(is.null(value) || is.list(value) || is(value, "List"))) {
        stop("'invalid 'value' in 'rowLoadings(<SingleCellExperiment>) <- value'")
    }

    for (i in seq_along(value)) {
        if (!identical(nrow(value[[i]]), nrow(x))) {
            stop("invalid 'value' in 'rowLoading(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "all elements in 'value' should have number of rows equal to 'nrow(x)'")
        }
        if (withDimnames && !is.null(rownames(value[[i]]))) {
            if (!identical(rownames(value[[i]]), rownames(x))) {
                warning("non-NULL 'rownames(value[[", i, "]])' should be the same as 'rownames(x)' for 'rowLoadings<-'.")
            }
        }
    }

    x <- .initialize_loadings(x)
    .set_internal_all(x, value,
                      getfun = int_elementMetadata,
                      setfun = `int_elementMetadata<-`,
                      key = .load_key,
                      convertfun = NULL,
                      xdimfun = nrow,
                      vdimfun = nrow,
                      funstr = "rowLoadings",
                      xdimstr = "nrow",
                      vdimstr = "rows")
})

#' @export
setGeneric("rowLoadingNames", function(x) standardGeneric("rowLoadingNames"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowLoadingNames", "SingleCellExperiment", function(x) {
    x <- .initialize_loadings(x)
    .get_internal_names(x, getfun = int_elementMetadata, key = .load_key)
})

#' @export
setGeneric("rowLoadingNames<-", function(x, value) standardGeneric("rowLoadingNames<-"))

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowLoadingNames", c("SingleCellExperiment", "character"), function(x, value) {
    x <- .initialize_loadings(x)
    .set_internal_names(x, value,
                        getfun = int_elementMetadata,
                        setfun = `int_elementMetadata<-`,
                        key = .load_key)
})

#' @export
setGeneric("rowLoading", function(x, type, withDimnames=TRUE) standardGeneric("rowLoading"))

#' @export
setMethod("rowLoading", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_loadings(x)
    .get_internal_missing(x,
                          basefun = rowLoading,
                          namefun = rowLoadingNames,
                          funstr = "rowLoading",
                          withDimnames = withDimnames)
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowLoading", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_loadings(x)
    out <- .get_internal_integer(x, type,
                                 getfun = int_elementMetadata,
                                 key = .load_key,
                                 funstr = "rowLoading",
                                 substr = "type")
    if (withDimnames) {
        rownames(out) <- rownames(x)
    }
    out
})

#' @export
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowLoading", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    x <- .initialize_loadings(x)
    out <- .get_internal_character(x, type,
                                   getfun = int_elementMetadata,
                                   key = .load_key,
                                   funstr = "rowLoading",
                                   substr = "type",
                                   namestr = "rowLoadingNames")
    if (withDimnames) {
        rownames(out) <- rownames(x)
    }
    out
})

#' @export
setGeneric("rowLoading<-", function(x, type, withDimnames=TRUE, ..., value) standardGeneric("rowLoading<-"))

#' @export
setReplaceMethod("rowLoading", c("SingleCellExperiment", "missing"), function(x, type, withDimnames=TRUE, ..., value) {
    x <- .initialize_loadings(x)
    .set_internal_missing(x, value,
                          withDimnames = withDimnames,
                          basefun = `rowLoading<-`,
                          namefun = rowLoadingNames)
})

#' @export
#' @importFrom S4Vectors isSingleNumber
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowLoading", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleNumber(type)) {
        stop("'type' must be a scalar in 'rowLoading(<SingleCellExperiment>, type=\"numeric\")' <- value'")
    }

    if (type > length(rowLoadingNames(x))) {
        stop("'type' out of bounds in 'rowLoading(<SingleCellExperiment>, type=\"numeric\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), nrow(x))) {
            stop("invalid 'value' in 'rowLoading(<SingleCellExperiment>, type=\"numeric\") <- value':\n  ",
                 "'value' should have number of rows equal to 'nrow(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), rownames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowLoading<-'.")
            }
        }
    }

    x <- .initialize_loadings(x)
    .set_internal_numeric(x, type, value,
                          getfun = int_elementMetadata,
                          setfun = `int_elementMetadata<-`,
                          key = .load_key,
                          convertfun = NULL,
                          xdimfun = nrow,
                          vdimfun = nrow,
                          funstr = "rowLoading",
                          xdimstr = "nrow",
                          vdimstr = "rows",
                          substr = "type")
})

#' @export
#' @importFrom S4Vectors isSingleString
#' @importFrom SingleCellExperiment int_elementMetadata int_elementMetadata<-
setReplaceMethod("rowLoading", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE, ..., value) {
    if (!isSingleString(type)) {
        stop("'type' must be a string in 'rowLoading(<SingleCellExperiment>, type=\"character\") <- value'")
    }

    if (!is.null(value)) {
        if (!identical(nrow(value), nrow(x))) {
            stop("invalid 'value' in 'rowLoading(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                 "'value' should have number of rows equal to 'nrow(x)'")
        }
        # Optional validation when withDimnames=TRUE
        if (withDimnames && !is.null(rownames(value))) {
            if (!identical(rownames(value), rownames(x))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowLoading<-'.")
            }
        }
    }

    x <- .initialize_loadings(x)
    .set_internal_character(x, type, value,
                            getfun = int_elementMetadata,
                            setfun = `int_elementMetadata<-`,
                            key = .load_key,
                            convertfun = NULL,
                            xdimfun = nrow,
                            vdimfun = nrow,
                            funstr = "rowLoading",
                            xdimstr = "nrow",
                            vdimstr = "rows",
                            substr = "type")
})
