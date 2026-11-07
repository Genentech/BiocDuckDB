#' DuckDBMatrix scuttle methods
#'
#' @description
#' scuttle methods for \link[DuckDBArray:DuckDBMatrix-class]{DuckDBMatrix} objects.
#'
#' @section QC Metrics Methods:
#' The following QC metrics methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{perCellQCMetrics(x, subsets = NULL, percent.top = integer(0),
#'         threshold = 0, ..., flatten = TRUE)}:}{
#'     Compute per-cell quality control metrics for a count matrix.
#'     \describe{
#'       \item{\code{subsets}}{named list of vectors specifying feature subsets
#'         (e.g., mitochondrial genes). Each vector can be character (feature names),
#'         logical, or integer (indices). Returns sum, detected, and percent for
#'         each subset.}
#'       \item{\code{percent.top}}{integer vector specifying the sizes of the top
#'         sets of high-abundance genes (e.g., \code{c(50, 100, 200)}). For each
#'         size, computes the percentage of total counts assigned to the most
#'         highly expressed genes in each cell. Uses SQL window functions for
#'         efficient computation.}
#'       \item{\code{threshold}}{numeric scalar specifying the detection threshold}
#'       \item{\code{flatten}}{logical scalar indicating whether nested DataFrames
#'         in the output should be flattened. If \code{TRUE} (default), subset
#'         and percent.top statistics are returned as top-level columns
#'         (e.g., \code{subsets_Mito_sum}, \code{percent.top_50}).
#'         If \code{FALSE}, nested DataFrames/matrices are returned.}
#'     }
#'   }
#'   \item{\code{perFeatureQCMetrics(x, subsets = NULL, threshold = 0, ..., flatten = TRUE)}:}{
#'     Compute per-feature quality control metrics for a count matrix.
#'     \describe{
#'       \item{\code{subsets}}{named list of vectors specifying cell subsets
#'         (e.g., control wells). Each vector can be character (cell names),
#'         logical, or integer (indices). Returns mean, detected, and ratio for
#'         each subset.}
#'       \item{\code{threshold}}{numeric scalar specifying the detection threshold}
#'       \item{\code{flatten}}{logical scalar indicating whether nested DataFrames
#'         in the output should be flattened. If \code{TRUE} (default), subset
#'         statistics are returned as top-level columns (e.g., \code{subsets_SetA_mean}).
#'         If \code{FALSE}, a nested \code{subsets} DataFrame is returned.}
#'     }
#'   }
#' }
#'
#' @section Normalization Methods:
#' The following normalization methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{librarySizeFactors(x)}:}{
#'     Define per-cell size factors from the library sizes.
#'   }
#'   \item{\code{geometricSizeFactors(x, pseudo.count = 1)}:}{
#'     Define per-cell size factors from the geometric mean of counts per cell.
#'     \describe{
#'       \item{\code{pseudo.count}}{numeric scalar specifying the pseudo-count to add
#'         during log-transformation}
#'     }
#'   }
#'   \item{\code{normalizeCounts(x, size.factors = NULL, log = TRUE, transform = c("log", "none", "asinh"),
#'         pseudo.count = 1, center.size.factors = TRUE, subset.row = NULL, ...)}:}{
#'     Normalizes counts by dividing by size factors.
#'     \describe{
#'       \item{\code{size.factors}}{numeric vector of cell-specific size factors}
#'       \item{\code{log}}{logical scalar indicating whether normalized values
#'         should be log2-transformed}
#'       \item{\code{transform}}{string specifying the transformation to apply
#'         to the normalized expression values.}
#'       \item{\code{pseudo.count}}{numeric scalar specifying the pseudo-count
#'         to add when \code{transform = "log"}}
#'       \item{\code{center.size.factors}}{logical scalar indicating whether
#'         size factors should be centered at unity before being used}
#'       \item{\code{subset.row}}{vector specifying the subset of rows of x for
#'         which to return normalized values}
#'     }
#'   }
#'   \item{\code{calculateTPM(x, lengths = NULL, ...)}:}{
#'     Calculates transcripts per million using SQL-optimized row-wise sweep.
#'     \describe{
#'       \item{\code{lengths}}{Numeric vector of gene lengths.}
#'     }
#'   }
#'   \item{\code{calculateCPM(x, size.factors = NULL, subset.row = NULL, ...)}:}{
#'     Calculates counts per million using DelayedArray-safe sweep operations.
#'     \describe{
#'       \item{\code{size.factors}}{numeric vector containing size factors to adjust
#'         the library sizes. If NULL, library sizes are used directly.}
#'       \item{\code{subset.row}}{vector specifying the subset of rows of x for
#'         which to return normalized values}
#'     }
#'   }
#'   \item{\code{calculateAverage(x, size.factors = NULL, subset.row = NULL, BPPARAM = SerialParam(), ...)}:}{
#'     Calculates average counts per feature using DelayedArray-safe sweep operations.
#'     \describe{
#'       \item{\code{size.factors}}{numeric vector containing size factors. If NULL,
#'         library size factors are computed automatically.}
#'       \item{\code{subset.row}}{vector specifying the subset of rows of x for
#'         which to return averages}
#'       \item{\code{BPPARAM}}{BiocParallelParam object for parallel computation}
#'     }
#'   }
#' }
#'
#' The following normalization methods will dispatch to optimized
#' DuckDBMatrix implementations of \code{normalizeCounts},
#' \code{librarySizeFactors}, \code{calculateCPM}, \code{calculateAverage},
#' and \code{calculateTPM}.
#'
#'
#' @section Aggregation Methods:
#' The following aggregation methods have optimized DuckDBMatrix implementations
#' that use SQL GROUP BY operations for efficient pseudo-bulk analysis:
#' \describe{
#'   \item{\code{numDetectedAcrossFeatures(x, ids, subset.row = NULL, subset.col = NULL,
#'         average = FALSE, threshold = 0, ...)}:}{
#'     Count detected features across feature sets for each cell.
#'     \describe{
#'       \item{\code{ids}}{factor or list specifying feature groupings}
#'       \item{\code{average}}{logical indicating whether to compute the proportion}
#'       \item{\code{threshold}}{threshold for detection}
#'     }
#'   }
#'   \item{\code{sumCountsAcrossFeatures(x, ids, subset.row = NULL, subset.col = NULL,
#'         average = FALSE, ...)}:}{
#'     Sum counts across feature sets for each cell.
#'     \describe{
#'       \item{\code{ids}}{factor or list specifying feature groupings}
#'       \item{\code{average}}{logical indicating whether to compute the average}
#'     }
#'   }
#'   \item{\code{summarizeAssayByGroup(x, ids, subset.row = NULL, subset.col = NULL,
#'         statistics = c("mean", "sum", "num.detected", "prop.detected"),
#'         store.number = "ncells", threshold = 0, ...)}:}{
#'     From an assay matrix, compute summary statistics for groups of cells.
#'     \describe{
#'       \item{\code{ids}}{factor specifying the group for each cell}
#'       \item{\code{statistics}}{character vector of statistics to compute}
#'       \item{\code{threshold}}{threshold for detection}
#'     }
#'   }
#' }
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[DuckDBArray:DuckDBMatrix-class]{DuckDBMatrix-class}} for the main class
#'   \item \code{\link[DuckDBArray]{DuckDBArray-matrixStats}} for the matrixStats methods
#'   \item \code{\link[scuttle]{perCellQCMetrics}} for the scuttle generic
#'   \item \code{\link[scuttle]{perFeatureQCMetrics}} for the scuttle generic
#'   \item \code{\link[scuttle]{librarySizeFactors}} for the scuttle generic
#'   \item \code{\link[scuttle]{normalizeCounts}} for the scuttle generic
#'   \item \code{\link[scuttle]{calculateCPM}} for the scuttle generic
#'   \item \code{\link[scuttle]{calculateAverage}} for the scuttle generic
#'   \item \code{\link[scuttle]{calculateTPM}} for the scuttle generic
#'   \item \code{\link[scuttle]{numDetectedAcrossFeatures}} for the scuttle generic
#'   \item \code{\link[scuttle]{sumCountsAcrossFeatures}} for the scuttle generic
#'   \item \code{\link[scuttle]{summarizeAssayByGroup}} for the scuttle generic
#' }
#'
#' @return
#' These methods mirror their \pkg{scuttle} generics.
#' \code{perCellQCMetrics()} and \code{perFeatureQCMetrics()} return a
#' \code{DataFrame} of QC metrics. \code{librarySizeFactors()},
#' \code{geometricSizeFactors()}, and \code{calculateAverage()} return a numeric
#' vector. \code{normalizeCounts()}, \code{calculateCPM()}, and
#' \code{calculateTPM()} return a matrix-like object of transformed values, and
#' \code{numDetectedAcrossFeatures()} and \code{sumCountsAcrossFeatures()} a
#' matrix-like object of per-group summaries. \code{summarizeAssayByGroup()}
#' returns a
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#' of per-group assay summaries.
#'
#' @aliases librarySizeFactors,DuckDBMatrix-method
#' @aliases geometricSizeFactors,DuckDBMatrix-method
#' @aliases normalizeCounts,DuckDBMatrix-method
#' @aliases calculateCPM,DuckDBMatrix-method
#' @aliases calculateAverage,DuckDBMatrix-method
#' @aliases calculateTPM,DuckDBMatrix-method
#' @aliases perCellQCMetrics,DuckDBMatrix-method
#' @aliases perFeatureQCMetrics,DuckDBMatrix-method
#' @aliases numDetectedAcrossFeatures,DuckDBMatrix-method
#' @aliases sumCountsAcrossFeatures,DuckDBMatrix-method
#' @aliases summarizeAssayByGroup,DuckDBMatrix-method
#'
#' @keywords utilities methods
#'
#' @examples
#' mat <- matrix(rpois(500, 5), 25, 20)
#' dimnames(mat) <- list(paste0("G", 1:25), paste0("C", 1:20))
#' names(dimnames(mat)) <- c("index1", "index2")
#' tmp <- tempfile()
#' writeParquet(mat, tmp)
#' pqmat <- DuckDBMatrix(tmp, datacol = "value",
#'     keycols = lapply(dimnames(mat), function(x) setNames(seq_along(x), x)))
#' head(librarySizeFactors(pqmat))
#' unlink(tmp, recursive = TRUE)
#'
#' @name DuckDBMatrix-scuttle
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### perCellQCMetrics
###

# Helper to convert subset specification to integer indices
# Follows scuttle::.subset2index pattern
.subset2index_DuckDBMatrix <- function(subset, target, byrow = TRUE)
{
    if (is.factor(subset)) {
        subset <- as.character(subset)
    }
    if (byrow) {
        dummy <- seq_len(nrow(target))
        names(dummy) <- rownames(target)
    } else {
        dummy <- seq_len(ncol(target))
        names(dummy) <- colnames(target)
    }

    if (!is.null(subset)) {
        subset <- dummy[subset]
        if (any(is.na(subset))) {
            stop("invalid subset indices specified")
        }
    } else {
        subset <- dummy
    }
    unname(subset)
}

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom dplyr collect group_by summarize
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom scuttle perCellQCMetrics
setMethod("perCellQCMetrics", "DuckDBMatrix",
function(x, subsets = NULL, percent.top = integer(0), threshold = 0, ...,
         flatten = TRUE, percent_top = NULL, detection_limit = NULL) {
    if (!is.null(percent_top)) {
        percent.top <- percent_top
    }
    if (!is.null(detection_limit)) {
        threshold <- detection_limit
    }
    # Validate subsets
    if (length(subsets) && is.null(names(subsets))) {
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index_DuckDBMatrix,
                      target = x, byrow = TRUE)

    # Validate and sort percent.top
    percent.top <- sort(as.integer(percent.top))

    # Use optimized single-query approach
    .perCellQCMetrics_DuckDBMatrix(x, subsets, percent.top, threshold, flatten)
})

# Single-query implementation for perCellQCMetrics using JOIN approach
#' @importFrom dplyr collect group_by inner_join mutate summarize n
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
.perCellQCMetrics_DuckDBMatrix <- function(x, subsets, percent.top, threshold, flatten)
{
    ncells <- ncol(x)
    ngenes <- nrow(x)
    cnames <- colnames(x)

    table <- x@seed@table
    row_key <- keynames(table)[1L]
    col_key <- keynames(table)[2L]
    row_keycol <- table@keycols[[1L]]
    col_keycol <- table@keycols[[2L]]
    datacol <- table@datacols[[1L]]
    fill <- x@seed@fill

    conn <- tblconn(table, select = FALSE)
    col_sym <- as.name(col_key)
    row_keyvals <- unname(row_keycol)
    col_keyvals <- unname(col_keycol)

    # Create subset membership dataframe for JOIN
    # Each row has: row_key, subset_name (or NA if not in any subset)
    if (length(subsets)) {
        subset_df_list <- lapply(names(subsets), function(nm) {
            data.frame(
                row_idx = row_keyvals[subsets[[nm]]],
                subset_name = nm,
                stringsAsFactors = FALSE
            )
        })
        subset_df <- do.call(rbind, subset_df_list)
        names(subset_df)[1L] <- row_key
    }

    # Total statistics query (always needed)
    total_aggr <- list(
        total_sum = call("sum", datacol, na.rm = TRUE),
        stored_n = call("n")
    )
    if (threshold == 0) {
        total_aggr$total_detected <- call("sum",
                                          call("as.integer", call("!=", datacol, 0L)),
                                          na.rm = TRUE)
    } else {
        total_aggr$total_detected <- call("sum",
                                          call("as.integer", call(">", datacol, threshold)),
                                          na.rm = TRUE)
    }

    total_result <- conn |>
        group_by(!!col_sym) |>
        summarize(!!!total_aggr, .groups = "drop") |>
        collect()

    # Convert integer64 to integer
    total_result$stored_n <- as.integer(total_result$stored_n)
    total_result$total_detected <- as.integer(total_result$total_detected)

    # Match results to output order
    idx <- match(col_keyvals, total_result[[col_key]])

    # Extract total statistics
    total_sum <- total_result$total_sum[idx]
    total_detected <- total_result$total_detected[idx]
    stored_n <- total_result$stored_n[idx]

    # Handle fill value for missing entries (sparse matrix)
    if (fill != 0) {
        missing_n <- ngenes - stored_n
        total_sum <- total_sum + fill * missing_n
        if (threshold < fill) {
            total_detected <- total_detected + missing_n
        }
    }

    # Handle cells with no entries in database
    total_sum[is.na(total_sum)] <- fill * ngenes
    total_detected[is.na(total_detected)] <- if (threshold < fill) ngenes else 0L

    full.info <- DataFrame(sum = total_sum,
                           detected = as.integer(total_detected),
                           row.names = cnames)

    # Compute percent.top if requested
    if (length(percent.top) > 0L) {
        pct_top <- .compute_percent_top_DuckDBMatrix(x, percent.top, total_sum)
        if (flatten) {
            # Flatten: add as separate columns
            for (i in seq_along(percent.top)) {
                col_name <- paste0("percent.top_", percent.top[i])
                full.info[[col_name]] <- pct_top[, i]
            }
        } else {
            full.info$percent.top <- pct_top
        }
    } else if (!flatten) {
        # Add empty percent.top matrix for compatibility with scuttle
        pct_top <- matrix(numeric(0), nrow = ncells, ncol = 0)
        dimnames(pct_top) <- list(NULL, NULL)
        full.info$percent.top <- pct_top
    }

    # Subset statistics using JOIN approach
    if (length(subsets)) {
        subset_aggr <- list(
            sub_sum = call("sum", datacol, na.rm = TRUE),
            stored_n = call("n")
        )
        if (threshold == 0) {
            subset_aggr$sub_detected <- call("sum",
                                             call("as.integer", call("!=", datacol, 0L)),
                                             na.rm = TRUE)
        } else {
            subset_aggr$sub_detected <- call("sum",
                                             call("as.integer", call(">", datacol, threshold)),
                                             na.rm = TRUE)
        }

        # subset_df's row_key column already contains exactly the rows that
        # belong to some subset, so an inner join on it is equivalent to (and
        # replaces) filter(row_key %in% that same set) |> left_join(...): no
        # row can pass that filter without also matching subset_df. This also
        # avoids ever materializing a literal IN (...) list of row keys in
        # the generated SQL, which does not scale to this package's
        # million-cell benchmarks.
        subset_result <- conn |>
            inner_join(subset_df, by = row_key, copy = TRUE) |>
            group_by(!!col_sym, !!as.name("subset_name")) |>
            summarize(!!!subset_aggr, .groups = "drop") |>
            collect()

        # Convert integer64 to integer
        subset_result$stored_n <- as.integer(subset_result$stored_n)
        subset_result$sub_detected <- as.integer(subset_result$sub_detected)

        # Pivot results into per-subset columns
        sub.info <- make_zero_col_DFrame(ncells)
        for (nm in names(subsets)) {
            subset_size <- length(subsets[[nm]])

            # Extract data for this subset
            sub_data <- subset_result[subset_result$subset_name == nm, , drop = FALSE]
            sub_idx <- match(col_keyvals, sub_data[[col_key]])

            sub_sum <- sub_data$sub_sum[sub_idx]
            sub_detected <- sub_data$sub_detected[sub_idx]

            # Handle fill value
            if (fill != 0) {
                sub_stored_n <- sub_data$stored_n[sub_idx]
                sub_stored_n[is.na(sub_stored_n)] <- 0L
                missing_sub_n <- subset_size - sub_stored_n
                sub_sum <- sub_sum + fill * missing_sub_n
                if (threshold < fill) {
                    sub_detected <- sub_detected + missing_sub_n
                }
            }

            # Handle NA values (cells with no entries for this subset)
            sub_sum[is.na(sub_sum)] <- fill * subset_size
            sub_detected[is.na(sub_detected)] <- if (threshold < fill) subset_size else 0L

            sub_percent <- sub_sum / total_sum * 100

            sub.info[[nm]] <- DataFrame(sum = sub_sum,
                                        detected = as.integer(sub_detected),
                                        percent = sub_percent)
        }
        full.info$subsets <- sub.info
        if (flatten) {
            full.info <- .flatten_nested_DuckDBMatrix(full.info)
        }
    } else if (!flatten) {
        # Add empty subsets DataFrame for compatibility with scuttle
        full.info$subsets <- make_zero_col_DFrame(ncells)
    }

    full.info
}

# Helper to flatten nested DataFrames (follows scuttle pattern)
#' @importFrom S4Vectors DataFrame
.flatten_nested_DuckDBMatrix <- function(x, name = "")
{
    if (!is.null(dim(x))) {
        if (name != "") {
            name <- paste0(name, "_")
        }
        names <- sprintf("%s%s", name, colnames(x))
        rn <- rownames(x)

        df <- vector("list", ncol(x))
        for (i in seq_along(df)) {
            df[[i]] <- .flatten_nested_DuckDBMatrix(x[, i], names[i])
        }
        if (length(df) > 0) {
            df <- do.call(cbind, df)
        } else {
            df <- DataFrame(x[, 0])
        }

        rownames(df) <- rn
    } else {
        df <- DataFrame(x)
        colnames(df) <- name
    }
    df
}

# Compute percent.top using SQL window functions
#' @importFrom dplyr collect filter group_by mutate sql summarize
#' @importFrom DuckDBDataFrame keynames tblconn
.compute_percent_top_DuckDBMatrix <- function(x, percent.top, total_sum)
{
    ncells <- ncol(x)

    table <- x@seed@table
    col_key <- keynames(table)[2L]
    col_keycol <- table@keycols[[2L]]
    datacol <- table@datacols[[1L]]
    datacol_name <- as.character(datacol)

    conn <- tblconn(table, select = FALSE)
    col_sym <- as.name(col_key)
    col_keyvals <- unname(col_keycol)
    max_top <- max(percent.top)

    # Use window function to rank genes by value within each cell
    # Then sum values for genes with rank <= threshold
    rank_col <- as.name("gene_rank")

    # Build aggregation expressions for each percent.top threshold
    aggr <- list()
    for (top_n in percent.top) {
        col_name <- paste0("top_", top_n)
        # Sum values where rank <= top_n
        aggr[[col_name]] <- call("sum",
                                 call("if_else", call("<=", rank_col, top_n), datacol, 0),
                                 na.rm = TRUE)
    }

    # SQL query with window function for ranking
    # ROW_NUMBER() OVER (PARTITION BY col_key ORDER BY value DESC)
    ranked_conn <- conn |>
        mutate(gene_rank = sql(paste0(
            "ROW_NUMBER() OVER (PARTITION BY \"", col_key,
            "\" ORDER BY \"", datacol_name, "\" DESC)"
        )))

    # Filter to only top N genes and aggregate
    result <- ranked_conn |>
        filter(!!rank_col <= max_top) |>
        group_by(!!col_sym) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    # Match results to output order
    idx <- match(col_keyvals, result[[col_key]])

    # Build output matrix
    pct_top <- matrix(NA_real_, nrow = ncells, ncol = length(percent.top))
    colnames(pct_top) <- as.character(percent.top)
    dimnames(pct_top) <- list(NULL, as.character(percent.top))

    for (i in seq_along(percent.top)) {
        col_name <- paste0("top_", percent.top[i])
        top_sum <- result[[col_name]][idx]
        top_sum[is.na(top_sum)] <- 0
        pct_top[, i] <- top_sum / total_sum * 100
    }

    pct_top
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### perFeatureQCMetrics
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom dplyr collect group_by summarize
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom scuttle perFeatureQCMetrics
setMethod("perFeatureQCMetrics", "DuckDBMatrix",
function(x, subsets = NULL, threshold = 0, ..., flatten = TRUE,
         detection_limit = NULL) {
    if (!is.null(detection_limit)) {
        threshold <- detection_limit
    }
    # Validate subsets
    if (length(subsets) && is.null(names(subsets))) {
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index_DuckDBMatrix,
                      target = x, byrow = FALSE)

    # Use optimized single-query approach
    .perFeatureQCMetrics_DuckDBMatrix(x, subsets, threshold, flatten)
})

# Single-query implementation for perFeatureQCMetrics using JOIN approach
#' @importFrom dplyr collect group_by inner_join summarize n
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
.perFeatureQCMetrics_DuckDBMatrix <- function(x, subsets, threshold, flatten)
{
    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    table <- x@seed@table
    row_key <- keynames(table)[1L]
    col_key <- keynames(table)[2L]
    row_keycol <- table@keycols[[1L]]
    col_keycol <- table@keycols[[2L]]
    datacol <- table@datacols[[1L]]
    fill <- x@seed@fill

    conn <- tblconn(table, select = FALSE)
    row_sym <- as.name(row_key)
    row_keyvals <- unname(row_keycol)
    col_keyvals <- unname(col_keycol)

    # Create subset membership dataframe for JOIN (by column)
    if (length(subsets)) {
        subset_df_list <- lapply(names(subsets), function(nm) {
            data.frame(
                col_idx = col_keyvals[subsets[[nm]]],
                subset_name = nm,
                stringsAsFactors = FALSE
            )
        })
        subset_df <- do.call(rbind, subset_df_list)
        names(subset_df)[1L] <- col_key
    }

    # Total statistics query (always needed)
    total_aggr <- list(
        total_sum = call("sum", datacol, na.rm = TRUE),
        stored_n = call("n")
    )
    if (threshold == 0) {
        total_aggr$total_detected <- call("sum",
                                          call("as.integer", call("!=", datacol, 0L)),
                                          na.rm = TRUE)
    } else {
        total_aggr$total_detected <- call("sum",
                                          call("as.integer", call(">", datacol, threshold)),
                                          na.rm = TRUE)
    }

    total_result <- conn |>
        group_by(!!row_sym) |>
        summarize(!!!total_aggr, .groups = "drop") |>
        collect()

    # Convert integer64 to integer
    total_result$stored_n <- as.integer(total_result$stored_n)
    total_result$total_detected <- as.integer(total_result$total_detected)

    # Match results to output order
    idx <- match(row_keyvals, total_result[[row_key]])

    # Extract total statistics
    total_sum <- total_result$total_sum[idx]
    total_detected <- total_result$total_detected[idx]
    stored_n <- total_result$stored_n[idx]

    # Handle fill value for missing entries (sparse matrix)
    if (fill != 0) {
        missing_n <- ncells - stored_n
        total_sum <- total_sum + fill * missing_n
        if (threshold < fill) {
            total_detected <- total_detected + missing_n
        }
    }

    # Handle genes with no entries in database
    total_sum[is.na(total_sum)] <- fill * ncells
    total_detected[is.na(total_detected)] <- if (threshold < fill) ncells else 0L

    total_mean <- total_sum / ncells
    detected_pct <- total_detected / ncells * 100

    full.info <- DataFrame(mean = total_mean,
                           detected = detected_pct,
                           row.names = rnames)

    # Subset statistics using JOIN approach
    if (length(subsets)) {
        subset_aggr <- list(
            sub_sum = call("sum", datacol, na.rm = TRUE),
            stored_n = call("n")
        )
        if (threshold == 0) {
            subset_aggr$sub_detected <- call("sum",
                                             call("as.integer", call("!=", datacol, 0L)),
                                             na.rm = TRUE)
        } else {
            subset_aggr$sub_detected <- call("sum",
                                             call("as.integer", call(">", datacol, threshold)),
                                             na.rm = TRUE)
        }

        # subset_df's col_key column already contains exactly the columns
        # that belong to some subset, so an inner join on it is equivalent to
        # (and replaces) filter(col_key %in% that same set) |> left_join(...)
        # (see the analogous comment in .perCellQCMetrics_DuckDBMatrix()).
        subset_result <- conn |>
            inner_join(subset_df, by = col_key, copy = TRUE) |>
            group_by(!!row_sym, !!as.name("subset_name")) |>
            summarize(!!!subset_aggr, .groups = "drop") |>
            collect()

        # Convert integer64 to integer
        subset_result$stored_n <- as.integer(subset_result$stored_n)
        subset_result$sub_detected <- as.integer(subset_result$sub_detected)

        # Pivot results into per-subset columns
        sub.info <- make_zero_col_DFrame(ngenes)
        for (nm in names(subsets)) {
            subset_size <- length(subsets[[nm]])

            # Extract data for this subset
            sub_data <- subset_result[subset_result$subset_name == nm, , drop = FALSE]
            sub_idx <- match(row_keyvals, sub_data[[row_key]])

            sub_sum <- sub_data$sub_sum[sub_idx]
            sub_detected <- sub_data$sub_detected[sub_idx]

            # Handle fill value
            if (fill != 0) {
                sub_stored_n <- sub_data$stored_n[sub_idx]
                sub_stored_n[is.na(sub_stored_n)] <- 0L
                missing_sub_n <- subset_size - sub_stored_n
                sub_sum <- sub_sum + fill * missing_sub_n
                if (threshold < fill) {
                    sub_detected <- sub_detected + missing_sub_n
                }
            }

            # Handle NA values (genes with no entries for this subset)
            sub_sum[is.na(sub_sum)] <- fill * subset_size
            sub_detected[is.na(sub_detected)] <- if (threshold < fill) subset_size else 0L

            sub_mean <- sub_sum / subset_size
            sub_detected_pct <- sub_detected / subset_size * 100
            sub_ratio <- sub_mean / total_mean

            sub.info[[nm]] <- DataFrame(mean = sub_mean,
                                        detected = sub_detected_pct,
                                        ratio = sub_ratio)
        }
        full.info$subsets <- sub.info
        if (flatten) {
            full.info <- .flatten_nested_DuckDBMatrix(full.info)
        }
    }

    full.info
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### librarySizeFactors
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom MatrixGenerics colSums
#' @importFrom scuttle librarySizeFactors
#' @importFrom stats setNames
setMethod("librarySizeFactors", "DuckDBMatrix",
function(x, subset.row = NULL, ..., subset_row = NULL) {
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    lib.sizes <- colSums(x)
    lib.sizes <- setNames(as.vector(lib.sizes), dimnames(lib.sizes)[[1L]])
    lib.sizes / mean(lib.sizes)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### geometricSizeFactors
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom MatrixGenerics colMeans
#' @importFrom scuttle geometricSizeFactors
#' @importFrom stats setNames
setMethod("geometricSizeFactors", "DuckDBMatrix",
function(x, subset.row = NULL, pseudo.count = 1, ...) {
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    means <- colMeans(log2(x + pseudo.count))
    means <- setNames(as.vector(means), dimnames(means)[[1L]])
    geo <- 2^means
    geo / mean(geo)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Fused library-size normalization (shared by normalizeCounts(),
### calculateAverage(), and calculateCPM() when size.factors is NULL)
###

# The common librarySizeFactors(x) -> colSums(x) -> sweep(x, 2L, ..., "/")
# chain costs a full round trip: colSums() forces a collect() of the
# per-cell sum into R, then sweep()'s left_join(conn, df, copy = TRUE)
# re-uploads that small vector as a DuckDB temp table and joins it back in --
# a second full table scan plus a temp-table upload, on top of the first
# scan colSums() already did.
#
# Replace that with ONE query using two chained window functions instead of
# a GROUP BY + R round trip + join:
#
#   size_factor_c = colsum_c / global_denom
#   normalized    = value / size_factor_c = value * global_denom / colsum_c
#
# where colsum_c = SUM(value) OVER (PARTITION BY cell) and global_denom is
# either a literal constant (calculateCPM()'s 1e6) or the mean library size
# across all cells (normalizeCounts()/calculateAverage()'s
# SUM(value) OVER () / ncells). This compiles to a single physical plan
# scanning the base table once, with the two window aggregates computed
# alongside it -- no join, no temp-table upload, no R-side collect() of the
# per-cell sums.
#
# Requires fill == 0, mirroring sweep()'s own hard requirement for the same
# architectural reason: a DuckDBMatrix's "fill" is a single scalar shared by
# every absent cell, but dividing by a per-cell size factor would turn a
# nonzero fill into a *different* normalized value in every column, which
# the seed's single-scalar fill slot cannot represent. This is not a new
# restriction introduced by the fusion: the two-pass chain it replaces
# already failed on a nonzero-fill DuckDBMatrix, because sweep() itself
# enforces the identical requirement (confirmed empirically:
# normalizeCounts() on a `fill != 0` DuckDBMatrix already errored with
# "must be a zero-filled array" before this change).
#
# Also reproduces the na.rm = FALSE NA-guard used elsewhere in this
# codebase (see .na_rm_guard() in DuckDBArray): librarySizeFactors()'s
# colSums() already makes any cell containing an NA have an NA library
# size, which poisons mean(lib.sizes) (na.rm = FALSE) and therefore every
# size factor and every output value in the WHOLE matrix, not just that
# cell's column (confirmed empirically against the pre-fusion
# implementation). A single global "any NA anywhere" flag reproduces this
# exactly, since once it's known no NA exists anywhere, no per-cell
# distinction is needed.
#' @importFrom DBI dbQuoteIdentifier
#' @importFrom dbplyr translate_sql
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom dplyr mutate sql
#' @importFrom stats setNames
.fused_normalize_by_libsize_DuckDBMatrix <-
function(x, mode = c("mean", "cpm"))
{
    mode <- match.arg(mode)
    if (x@seed@fill != 0) {
        stop("must be a zero-filled array")
    }
    table <- x@seed@table
    if (length(table@datacols) != 1L) {
        stop("requires a single datacol")
    }
    col_key <- keynames(table)[2L]
    ncells <- ncol(x)

    conn <- tblconn(table, select = FALSE)
    db_conn <- conn$src$con
    datacol <- table@datacols[[1L]]
    datacol_sql <- as.character(translate_sql(!!datacol, con = db_conn))
    col_key_sql <- as.character(dbQuoteIdentifier(db_conn, col_key))

    global_denom_sql <- if (mode == "cpm") {
        "1000000.0"
    } else {
        sprintf("(SUM(%s) OVER () / %s)", datacol_sql, sprintf("%.17g", ncells))
    }
    has_na_sql <- sprintf("BOOL_OR((%s) IS NULL) OVER ()", datacol_sql)
    colsum_cell_sql <- sprintf("SUM(%s) OVER (PARTITION BY %s)", datacol_sql, col_key_sql)
    normalized_sql <- sprintf(
        "CASE WHEN %s THEN NULL ELSE (%s) * (%s) / (%s) END",
        has_na_sql, datacol_sql, global_denom_sql, colsum_cell_sql
    )

    new_conn <- mutate(conn, normalized_value = sql(normalized_sql))
    new_datacols <- setNames(as.expression(as.name("normalized_value")),
                             names(table@datacols))
    new_table <- replaceSlots(table, conn = new_conn, datacols = new_datacols,
                              check = FALSE)
    new_seed <- replaceSlots(x@seed, table = new_table, check = FALSE)
    replaceSlots(x, seed = new_seed, check = FALSE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### normalizeCounts
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom MatrixGenerics colSums
#' @importFrom scuttle normalizeCounts
setMethod("normalizeCounts", "DuckDBMatrix",
function(x, size.factors = NULL, log = TRUE, transform = c("log", "none", "asinh"),
         pseudo.count = 1, center.size.factors = TRUE, subset.row = NULL, ...,
         size_factors = NULL, center_size_factors = NULL, pseudo_count = NULL,
         subset_row = NULL) {
    if (!is.null(size_factors)) {
        size.factors <- size_factors
    }
    if (!is.null(center_size_factors)) {
        center.size.factors <- center_size_factors
    }
    if (!is.null(pseudo_count)) {
        pseudo.count <- pseudo_count
    }
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    transform <- match.arg(transform)
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    if (is.null(size.factors)) {
        out <- .fused_normalize_by_libsize_DuckDBMatrix(x, mode = "mean")
    } else {
        if (center.size.factors) {
            size.factors <- size.factors / mean(size.factors)
        }
        out <- sweep(x, 2L, size.factors, FUN = "/")
    }
    if (transform == "log" && log) {
        out <- log(out + pseudo.count) / log(2)
    } else if (transform == "asinh") {
        out <- asinh(out) / log(2)
    }
    out
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### calculateCPM
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom MatrixGenerics colSums
#' @importFrom scuttle calculateCPM
setMethod("calculateCPM", "DuckDBMatrix",
function(x, size.factors = NULL, subset.row = NULL, ...,
         size_factors = NULL, subset_row = NULL) {
    if (!is.null(size_factors)) {
        size.factors <- size_factors
    }
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    if (is.null(size.factors)) {
        .fused_normalize_by_libsize_DuckDBMatrix(x, mode = "cpm")
    } else {
        lib.sizes <- colSums(x) / 1e6
        lib.sizes <- size.factors / mean(size.factors) * mean(lib.sizes)
        sweep(x, 2L, lib.sizes, FUN = "/")
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### calculateAverage
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom MatrixGenerics rowMeans
#' @importFrom BiocParallel SerialParam
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom scuttle calculateAverage
setMethod("calculateAverage", "DuckDBMatrix",
function(x, size.factors = NULL, subset.row = NULL, BPPARAM = SerialParam(), ...,
         size_factors = NULL, subset_row = NULL) {
    if (!is.null(size_factors)) {
        size.factors <- size_factors
    }
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    if (is.null(size.factors)) {
        x <- .fused_normalize_by_libsize_DuckDBMatrix(x, mode = "mean")
    } else {
        size.factors <- size.factors / mean(size.factors)
        x <- sweep(x, 2L, size.factors, FUN = "/")
    }
    oldbp <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldbp))

    rowMeans(x)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### calculateTPM
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom scuttle calculateTPM
setMethod("calculateTPM", "DuckDBMatrix",
function(x, lengths = NULL, size.factors = NULL, subset.row = NULL, ...,
         size_factors = NULL, subset_row = NULL) {
    if (!is.null(size_factors)) {
        size.factors <- size_factors
    }
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
        if (!is.null(lengths)) {
            lengths <- lengths[subset.row]
        }
    }
    if (!is.null(lengths)) {
        x <- sweep(x, 1L, lengths, FUN = "/")
    }
    calculateCPM(x, size.factors = size.factors, ...)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Group aggregation helper functions
###

.processIds <- function(x, ids, subset.col) {
    if (length(ids) != ncol(x)) {
        stop("length of 'ids' and 'ncol(x)' are not equal")
    }
    ids <- as.factor(ids)
    valid <- !is.na(ids)
    if (!is.null(subset.col)) {
        valid <- valid & (seq_along(ids) %in% subset.col)
    }
    list(ids = ids, valid = valid)
}

#' @importFrom DuckDBDataFrame keynames
.processFeatureIds <- function(x, ids, subset.row) {
    table <- x@seed@table
    row_key <- keynames(table)[1L]
    row_keycol <- table@keycols[[1L]]
    if (is.list(ids)) {
        feature_df <- do.call(rbind, lapply(names(ids), function(nm) {
            data.frame(row_idx = as.integer(ids[[nm]]), feature_group = nm,
                       stringsAsFactors = FALSE)
        }))
        feature_df[[row_key]] <- unname(row_keycol[feature_df$row_idx])
    } else {
        if (length(ids) != nrow(x)) {
            stop("'ids' should be of length equal to 'nrow(x)'")
        }
        ids <- as.factor(ids)
        valid <- !is.na(ids)
        feature_df <- data.frame(
            row_idx = which(valid),
            feature_group = as.character(ids[valid]),
            stringsAsFactors = FALSE
        )
        feature_df[[row_key]] <- unname(row_keycol[valid])
    }
    if (!is.null(subset.row)) {
        feature_df <- feature_df[feature_df$row_idx %in% subset.row, , drop = FALSE]
    }
    feature_df
}

.pivotToMatrix <- function(data, value_col, row_col, col_col, nrow, ncol,
                           rnames, cnames, default = 0) {
    mat <- matrix(default, nrow = nrow, ncol = ncol, dimnames = list(rnames, cnames))
    idx <- cbind(match(data[[row_col]], rnames), match(data[[col_col]], cnames))
    mat[idx] <- data[[value_col]]
    mat
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### numDetectedAcrossFeatures
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom dplyr collect group_by inner_join n summarize
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom scuttle numDetectedAcrossFeatures
setMethod("numDetectedAcrossFeatures", "DuckDBMatrix",
function(x, ids, subset.row = NULL, subset.col = NULL, average = FALSE,
         threshold = 0, ..., subset_row = NULL, subset_col = NULL,
         detection_limit = NULL) {
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset_col)) {
        subset.col <- subset_col
    }
    if (!is.null(detection_limit)) {
        threshold <- detection_limit
    }
    feature_df <- .processFeatureIds(x, ids, subset.row)

    if (!is.null(subset.col)) {
        x <- x[, subset.col, drop = FALSE]
    }

    feature_col <- "feature_group"
    groups <- unique(feature_df[[feature_col]])
    ngroups <- length(groups)
    ncells <- ncol(x)
    cnames <- colnames(x)
    group_sizes <- table(feature_df[[feature_col]])[groups]

    table <- x@seed@table
    row_key <- keynames(table)[1L]
    col_key <- keynames(table)[2L]
    col_keycol <- table@keycols[[2L]]
    col_keyvals <- unname(col_keycol)

    conn <- tblconn(table, select = FALSE)
    fill <- x@seed@fill
    datacol <- table@datacols[[1L]]

    detect_expr <- if (threshold == 0)
                       call("countif", call("!=", call("(", datacol), 0L))
                   else
                       call("countif", call(">", call("(", datacol), threshold))
    aggr <- list(detected = detect_expr, stored_n = call("n"))

    # feature_df's row_key column already contains exactly the rows assigned
    # to a feature group, so an inner join on it is equivalent to (and
    # replaces) filter(row_key %in% that same set) |> left_join(...) (see the
    # analogous comment in .perCellQCMetrics_DuckDBMatrix()).
    result <- conn |>
        inner_join(feature_df[, c(row_key, feature_col)], by = row_key, copy = TRUE) |>
        group_by(!!as.name(feature_col), !!as.name(col_key)) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    out <- .pivotToMatrix(result, "detected", feature_col, col_key,
                          ngroups, ncells, groups, col_keyvals, default = 0L)
    colnames(out) <- cnames

    if (fill != 0 && threshold < fill) {
        count_mat <- .pivotToMatrix(result, "stored_n", feature_col, col_key,
                                    ngroups, ncells, groups, col_keyvals)
        out <- out + (as.numeric(group_sizes) - count_mat)
    }

    if (average) {
        out <- out / as.numeric(group_sizes)
    }

    out
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sumCountsAcrossFeatures
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom dplyr collect group_by inner_join n summarize
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom scuttle sumCountsAcrossFeatures
setMethod("sumCountsAcrossFeatures", "DuckDBMatrix",
function(x, ids, subset.row = NULL, subset.col = NULL, average = FALSE, ...,
         subset_row = NULL, subset_col = NULL) {
    if (!is.null(subset_row)) {
        subset.row <- subset_row
    }
    if (!is.null(subset_col)) {
        subset.col <- subset_col
    }
    feature_df <- .processFeatureIds(x, ids, subset.row)

    if (!is.null(subset.col)) {
        x <- x[, subset.col, drop = FALSE]
    }

    feature_col <- "feature_group"
    groups <- unique(feature_df[[feature_col]])
    ngroups <- length(groups)
    ncells <- ncol(x)
    cnames <- colnames(x)
    group_sizes <- table(feature_df[[feature_col]])[groups]

    table <- x@seed@table
    row_key <- keynames(table)[1L]
    col_key <- keynames(table)[2L]
    col_keycol <- table@keycols[[2L]]
    col_keyvals <- unname(col_keycol)

    conn <- tblconn(table, select = FALSE)
    fill <- x@seed@fill
    datacol <- table@datacols[[1L]]

    aggr <- list(total = call("sum", datacol, na.rm = TRUE), stored_n = call("n"))

    # See the analogous comment in numDetectedAcrossFeatures(): feature_df's
    # row_key column is already exactly the filtered row set, so inner_join
    # replaces filter(...) |> left_join(...).
    result <- conn |>
        inner_join(feature_df[, c(row_key, feature_col)], by = row_key, copy = TRUE) |>
        group_by(!!as.name(feature_col), !!as.name(col_key)) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    out <- .pivotToMatrix(result, "total", feature_col, col_key,
                          ngroups, ncells, groups, col_keyvals)
    colnames(out) <- cnames

    if (fill != 0) {
        count_mat <- .pivotToMatrix(result, "stored_n", feature_col, col_key,
                                    ngroups, ncells, groups, col_keyvals)
        out <- out + (fill * (as.numeric(group_sizes) - count_mat))
    }

    if (average) {
        out <- out / as.numeric(group_sizes)
    }

    out
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarizeAssayByGroup
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom dplyr collect group_by inner_join n summarize
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle summarizeAssayByGroup
#' @importFrom SummarizedExperiment SummarizedExperiment
setMethod("summarizeAssayByGroup", "DuckDBMatrix",
function(x, ids, subset.row = NULL, subset.col = NULL,
         statistics = c("mean", "sum", "num.detected", "prop.detected"),
         store.number = "ncells", threshold = 0, ...) {
    statistics <- match.arg(statistics, several.ok = TRUE)

    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    id_info <- .processIds(x, ids, subset.col)
    ids <- id_info$ids
    valid <- id_info$valid

    groups <- levels(ids)[levels(ids) %in% ids[valid]]
    ngroups <- length(groups)
    ngenes <- nrow(x)
    rnames <- rownames(x)
    group_sizes <- table(ids[valid])[groups]

    table <- x@seed@table
    row_key <- keynames(table)[1L]
    col_key <- keynames(table)[2L]
    group_col <- "group_id"
    row_keycol <- table@keycols[[1L]]
    col_keycol <- table@keycols[[2L]]
    col_vals <- unname(col_keycol[valid])
    group_df <- data.frame(col_key = col_vals, group_id = as.character(ids[valid]),
                           stringsAsFactors = FALSE)
    names(group_df)[1L] <- col_key

    conn <- tblconn(table, select = FALSE)
    fill <- x@seed@fill
    datacol <- table@datacols[[1L]]

    aggr <- list(
        total = if ("sum" %in% statistics || "mean" %in% statistics)
                    call("sum", datacol, na.rm = TRUE) else NULL,
        detected = if ("num.detected" %in% statistics || "prop.detected" %in% statistics)
                       if (threshold == 0)
                           call("countif", call("!=", call("(", datacol), 0L))
                       else
                           call("countif", call(">", call("(", datacol), threshold))
                   else NULL,
        stored_n = call("n")
    )
    aggr <- aggr[!vapply(aggr, is.null, logical(1L))]

    # group_df's col_key column already contains exactly the (valid,
    # subset.col-filtered) columns being grouped, so an inner join on it is
    # equivalent to (and replaces) filter(col_key %in% that same set) |>
    # left_join(...) (see the analogous comment in
    # .perCellQCMetrics_DuckDBMatrix()).
    result <- conn |>
        inner_join(group_df, by = col_key, copy = TRUE) |>
        group_by(!!as.name(row_key), !!as.name(group_col)) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    row_keyvals <- unname(row_keycol)
    collected <- list()

    if ("sum" %in% statistics || "mean" %in% statistics) {
        stored_sum <- .pivotToMatrix(result, "total", row_key, "group_id",
                                     ngenes, ngroups, row_keyvals, groups)
        rownames(stored_sum) <- rnames
        if (fill != 0) {
            count_mat <- .pivotToMatrix(result, "stored_n", row_key, "group_id",
                                        ngenes, ngroups, row_keyvals, groups)
            stored_sum <- stored_sum +
                fill * t(as.numeric(group_sizes) - t(count_mat))
        }
        if ("sum" %in% statistics) {
            collected$sum <- stored_sum
        }
        if ("mean" %in% statistics) {
            collected$mean <- t(t(stored_sum) / as.numeric(group_sizes))
        }
    }

    if ("num.detected" %in% statistics || "prop.detected" %in% statistics) {
        stored_detected <- .pivotToMatrix(result, "detected", row_key, "group_id",
                                          ngenes, ngroups, row_keyvals, groups, default = 0L)
        rownames(stored_detected) <- rnames
        if (fill != 0 && threshold < fill) {
            count_mat <- .pivotToMatrix(result, "stored_n", row_key, "group_id",
                                        ngenes, ngroups, row_keyvals, groups)
            stored_detected <- stored_detected + t(as.numeric(group_sizes) - t(count_mat))
        }
        if ("num.detected" %in% statistics) {
            collected$num.detected <- stored_detected
        }
        if ("prop.detected" %in% statistics) {
            collected$prop.detected <- t(t(stored_detected) / as.numeric(group_sizes))
        }
    }

    coldata <- DataFrame(ids = groups, row.names = groups)
    if (!is.null(store.number)) {
        coldata[[store.number]] <- as.integer(group_sizes)
    }

    SummarizedExperiment(collected[statistics], colData = coldata)
})
