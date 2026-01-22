#' DuckDBMatrix scuttle methods
#'
#' @description
#' scuttle methods for \linkS4class{DuckDBMatrix} objects.
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
#' }
#'
#' The following normalization methods will dispatch to optimized
#' DuckDBMatrix implementations of \code{normalizeCounts},
#' \code{librarySizeFactors}, and \code{calculateTPM}:
#' \describe{
#'   \item{\code{calculateAverage(x, ...)}:}{
#'     Calculates average normalized expression per feature.
#'   }
#'   \item{\code{calculateCPM(x, ...)}:}{
#'     Calculates counts per million.
#'   }
#' }
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
#'   \item \code{\link{DuckDBMatrix-class}} for the main class
#'   \item \code{\link{DuckDBArray-matrixStats}} for the matrixStats methods
#'   \item \code{\link[scuttle]{perCellQCMetrics}} for the scuttle generic
#'   \item \code{\link[scuttle]{perFeatureQCMetrics}} for the scuttle generic
#'   \item \code{\link[scuttle]{librarySizeFactors}} for the scuttle generic
#'   \item \code{\link[scuttle]{normalizeCounts}} for the scuttle generic
#'   \item \code{\link[scuttle]{calculateTPM}} for the scuttle generic
#'   \item \code{\link[scuttle]{numDetectedAcrossFeatures}} for the scuttle generic
#'   \item \code{\link[scuttle]{sumCountsAcrossFeatures}} for the scuttle generic
#'   \item \code{\link[scuttle]{summarizeAssayByGroup}} for the scuttle generic
#' }
#'
#' @aliases
#' librarySizeFactors,DuckDBMatrix-method
#' geometricSizeFactors,DuckDBMatrix-method
#' normalizeCounts,DuckDBMatrix-method
#' calculateTPM,DuckDBMatrix-method
#' perCellQCMetrics,DuckDBMatrix-method
#' perFeatureQCMetrics,DuckDBMatrix-method
#' numDetectedAcrossFeatures,DuckDBMatrix-method
#' sumCountsAcrossFeatures,DuckDBMatrix-method
#' summarizeAssayByGroup,DuckDBMatrix-method
#'
#' @keywords utilities methods
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
         flatten = TRUE) {
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
#' @importFrom dplyr collect filter group_by left_join mutate summarize n
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

        # Filter to only rows in subsets, then join with subset membership
        all_subset_rows <- unique(unlist(lapply(subsets, function(s) row_keyvals[s])))

        subset_result <- conn |>
            filter(!!as.name(row_key) %in% local(all_subset_rows)) |>
            left_join(subset_df, by = row_key, copy = TRUE) |>
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
function(x, subsets = NULL, threshold = 0, ..., flatten = TRUE) {
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
#' @importFrom dplyr collect filter group_by left_join summarize n
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

        # Filter to only columns in subsets, then join with subset membership
        all_subset_cols <- unique(unlist(lapply(subsets, function(s) col_keyvals[s])))

        subset_result <- conn |>
            filter(!!as.name(col_key) %in% local(all_subset_cols)) |>
            left_join(subset_df, by = col_key, copy = TRUE) |>
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
function(x, subset.row = NULL, ...) {
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
### normalizeCounts
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom MatrixGenerics colSums
#' @importFrom scuttle normalizeCounts
setMethod("normalizeCounts", "DuckDBMatrix",
function(x, size.factors = NULL, log = TRUE, transform = c("log", "none", "asinh"),
         pseudo.count = 1, center.size.factors = TRUE, subset.row = NULL, ...) {
    transform <- match.arg(transform)
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x)
    } else if (center.size.factors) {
        size.factors <- size.factors / mean(size.factors)
    }
    out <- sweep(x, 2L, size.factors, FUN = "/")
    if (transform == "log" && log) {
        out <- log(out + pseudo.count) / log(2)
    } else if (transform == "asinh") {
        out <- asinh(out) / log(2)
    }
    out
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### calculateTPM
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom DelayedArray sweep
#' @importFrom scuttle calculateCPM calculateTPM
setMethod("calculateTPM", "DuckDBMatrix",
function(x, lengths = NULL, size.factors = NULL, subset.row = NULL, ...) {
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
#' @importFrom dplyr collect filter group_by left_join n summarize
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom scuttle numDetectedAcrossFeatures
setMethod("numDetectedAcrossFeatures", "DuckDBMatrix",
function(x, ids, subset.row = NULL, subset.col = NULL, average = FALSE,
         threshold = 0, ...) {
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

    result <- conn |>
        filter(!!as.name(row_key) %in% local(feature_df[[row_key]])) |>
        left_join(feature_df[, c(row_key, feature_col)], by = row_key, copy = TRUE) |>
        group_by(!!as.name(feature_col), !!as.name(col_key)) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    out <- .pivotToMatrix(result, "detected", feature_col, col_key,
                          ngroups, ncells, groups, col_keyvals, default = 0L)
    colnames(out) <- cnames

    if (fill != 0 && threshold < fill) {
        count_mat <- .pivotToMatrix(result, "stored_n", feature_col, col_key,
                                    ngroups, ncells, groups, col_keyvals)
        for (i in seq_len(ngroups)) {
            out[i, ] <- out[i, ] + (group_sizes[i] - count_mat[i, ])
        }
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
#' @importFrom dplyr collect filter group_by left_join n summarize
#' @importFrom DuckDBDataFrame keynames tblconn
#' @importFrom scuttle sumCountsAcrossFeatures
setMethod("sumCountsAcrossFeatures", "DuckDBMatrix",
function(x, ids, subset.row = NULL, subset.col = NULL, average = FALSE, ...) {
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

    result <- conn |>
        filter(!!as.name(row_key) %in% local(feature_df[[row_key]])) |>
        left_join(feature_df[, c(row_key, feature_col)], by = row_key, copy = TRUE) |>
        group_by(!!as.name(feature_col), !!as.name(col_key)) |>
        summarize(!!!aggr, .groups = "drop") |>
        collect()

    out <- .pivotToMatrix(result, "total", feature_col, col_key,
                          ngroups, ncells, groups, col_keyvals)
    colnames(out) <- cnames

    if (fill != 0) {
        count_mat <- .pivotToMatrix(result, "stored_n", feature_col, col_key,
                                    ngroups, ncells, groups, col_keyvals)
        for (i in seq_len(ngroups)) {
            out[i, ] <- out[i, ] + (fill * (group_sizes[i] - count_mat[i, ]))
        }
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
#' @importFrom dplyr collect filter group_by left_join n summarize
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

    result <- conn |>
        filter(!!as.name(col_key) %in% local(group_df[[col_key]])) |>
        left_join(group_df, by = col_key, copy = TRUE) |>
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
            for (j in seq_len(ngroups)) {
                stored_sum[, j] <- stored_sum[, j] + (fill * (group_sizes[j] - count_mat[, j]))
            }
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
            for (j in seq_len(ngroups)) {
                stored_detected[, j] <- stored_detected[, j] + (group_sizes[j] - count_mat[, j])
            }
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
