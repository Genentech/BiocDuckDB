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
#' @aliases
#' rowLoading
#' rowLoading,SingleCellExperiment,numeric-method
#' rowLoading,SingleCellExperiment,character-method
#' rowLoadings
#' rowLoadings,SingleCellExperiment-method
#' rowLoadingNames
#' rowLoadingNames,SingleCellExperiment-method
#' rowLoading<-
#' rowLoading<-,SingleCellExperiment,numeric-method
#' rowLoading<-,SingleCellExperiment,character-method
#' rowLoadings<-
#' rowLoadings<-,SingleCellExperiment-method
#' rowLoadingNames<-
#' rowLoadingNames<-,SingleCellExperiment,character-method
#'
#' @name SingleCellExperiment-loadings
NULL

#' @export
setGeneric("rowLoadings", function(x, withDimnames=TRUE) standardGeneric("rowLoadings"))

#' @export
#' @importClassesFrom S4Vectors SimpleList
#' @importFrom S4Vectors endoapply make_zero_col_DFrame
#' @importFrom SingleCellExperiment int_elementMetadata
setMethod("rowLoadings", "SingleCellExperiment", function(x, withDimnames=TRUE) {
    loadings <- int_elementMetadata(x)[["rowLoadings"]]
    if (is.null(loadings)) {
        loadings <- make_zero_col_DFrame(nrow(x))
    }
    loadings <- as(loadings, "SimpleList")
    if (withDimnames) {
        loadings <- endoapply(loadings, function(mat) {
            rownames(mat) <- rownames(x)
            mat
        })
    }
    loadings
})

#' @export
setGeneric("rowLoadings<-", function(x, withDimnames=TRUE, ..., value) standardGeneric("rowLoadings<-"))

#' @export
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors DataFrame I mcols mcols<- metadata metadata<-
setReplaceMethod("rowLoadings", "SingleCellExperiment", function(x, withDimnames=TRUE, ..., value) {
    if (!(is.null(value) || is.list(value) || is(value, "List"))) {
        stop("'invalid 'value' in 'rowLoadings(<SingleCellExperiment>) <- value'")
    }
    if (length(value) == 0) {
        loadings <- NULL
    } else {
        for (i in seq_len(length(value))) {
            if (!identical(nrow(value[[i]]), nrow(x))) {
                stop("invalid 'value' in 'rowLoading(<SingleCellExperiment>, type=\"character\") <- value':\n  ",
                     "all elements in 'value' should have number of rows equal to 'nrow(x)'")
            }
            if (withDimnames && !is.null(rownames(value[[i]]))) {
                if (!identical(rownames(x), rownames(value[[i]]))) {
                    warning("non-NULL 'rownames(value[[", i, "]])' should be the same as 'rownames(x)' for 'rowLoadings<-'.")
                }
            }
        }

        if (is.null(names(value))) {
            names(value) <- paste0("loadings", seq_along(value))
        }
        loadings <- do.call(DataFrame, c(lapply(value, I), list(row.names=NULL, check.names=FALSE)))

        if (is(value, "Annotated")) {
            metadata(loadings) <- metadata(value)
        }

        if (is(value, "Vector")) {
            mcols(loadings) <- mcols(value)
        }
    }
    int_elementMetadata(x)[["rowLoadings"]] <- loadings
    x
})

#' @export
setGeneric("rowLoadingNames", function(x) standardGeneric("rowLoadingNames"))

#' @export
setMethod("rowLoadingNames", "SingleCellExperiment", function(x) names(rowLoadings(x)))

#' @export
setGeneric("rowLoadingNames<-", function(x, value) standardGeneric("rowLoadingNames<-"))

#' @export
setReplaceMethod("rowLoadingNames", c("SingleCellExperiment", "character"), function(x, value) {
    if (length(value) != length(rowLoadings(x))) {
        stop("invalid 'value' in 'rowLoadingNames(<SingleCellExperiment>) <- value':\n  ",
             "'value' should have length equal to the number of results in 'x'")
    }
    loadings <- rowLoadings(x)
    names(loadings) <- value
    rowLoadings(x) <- loadings
    x
})

#' @export
setGeneric("rowLoading", function(x, type, withDimnames=TRUE) standardGeneric("rowLoading"))

#' @export
#' @importFrom S4Vectors isSingleNumber
setMethod("rowLoading", c("SingleCellExperiment", "numeric"), function(x, type, withDimnames=TRUE) {
    if (!isSingleNumber(type)) {
        stop("'type' must be a scalar in 'rowLoading(<SingleCellExperiment>, type=\"numeric\")'")
    }
    if (type > length(rowLoadingNames(x))) {
        stop("invalid subscript 'type' in 'rowLoading(<SingleCellExperiment>, type=\"numeric\")'")
    }
    rowLoadings(x, withDimnames=withDimnames)[[type]]
})

#' @export
#' @importFrom S4Vectors isSingleString
setMethod("rowLoading", c("SingleCellExperiment", "character"), function(x, type, withDimnames=TRUE) {
    if (!isSingleString(type)) {
        stop("'type' must be a string in 'rowLoading(<SingleCellExperiment>, type=\"character\")'")
    }
    if (!type %in% rowLoadingNames(x)) {
        stop("invalid subscript 'type' in 'rowLoadingNames(<SingleCellExperiment>, type=\"character\")'")
    }
    rowLoadings(x, withDimnames=withDimnames)[[type]]
})

#' @export
setGeneric("rowLoading<-", function(x, type, withDimnames=TRUE, ..., value) standardGeneric("rowLoading<-"))

#' @export
#' @importFrom S4Vectors isSingleNumber
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
            if (!identical(rownames(x), rownames(value))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowLoading<-'.")
            }
        }
    }
    loadings <- rowLoadings(x, withDimnames=FALSE)
    loadings[[type]] <- value
    rowLoadings(x, withDimnames=FALSE) <- loadings
    x
})

#' @export
#' @importFrom S4Vectors isSingleString
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
            if (!identical(rownames(x), rownames(value))) {
                warning("non-NULL 'rownames(value)' should be the same as 'rownames(x)' for 'rowLoading<-'.")
            }
        }
    }
    loadings <- rowLoadings(x, withDimnames=FALSE)
    loadings[[type]] <- value
    rowLoadings(x, withDimnames=FALSE) <- loadings
    x
})
