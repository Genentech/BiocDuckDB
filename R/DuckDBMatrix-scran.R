#' DuckDBMatrix scran methods
#'
#' @description
#' scran methods for \linkS4class{DuckDBMatrix} objects.
#'
#' @section Pairwise Correlation Methods:
#' The following pairwise correlation methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{correlatePairs(x, subset.row = NULL, pairings = NULL,
#'     use.names = TRUE, BPPARAM = SerialParam())}:}{
#'     Compute pairwise Pearson correlations between genes using sparse-aware
#'     SQL aggregation. Only \code{fill = 0} is supported.
#'     \describe{
#'       \item{\code{subset.row}}{rows (genes) to correlate; required for >1000 genes}
#'       \item{\code{pairings}}{optional matrix specifying specific gene pairs}
#'       \item{\code{use.names}}{whether to use row names in output}
#'     }
#'   }
#' }
#'
#' @section Variance Modelling Methods:
#' The following variance modelling methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{modelGeneVar(x, block = NULL, design = NULL, subset.row = NULL,
#'     subset.fit = NULL, ..., equiweight = TRUE, method = "fisher",
#'     BPPARAM = SerialParam())}:}{
#'     Model the variance of log-expression profiles for each gene,
#'     decomposing it into technical and biological components.
#'     Uses SQL aggregation for computing per-gene means and variances.
#'     \describe{
#'       \item{\code{block}}{factor specifying blocking levels for each cell}
#'       \item{\code{design}}{design matrix for blocking (not supported with block)}
#'       \item{\code{subset.row}}{rows for which to model the variance}
#'       \item{\code{subset.fit}}{rows to be used for trend fitting}
#'       \item{\code{equiweight}}{whether blocks should be weighted equally}
#'       \item{\code{method}}{method for combining p-values across blocks}
#'     }
#'   }
#'   \item{\code{modelGeneVarByPoisson(x, size.factors = NULL, block = NULL,
#'     design = NULL, subset.row = NULL, ..., equiweight = TRUE, method = "fisher",
#'     BPPARAM = SerialParam())}:}{
#'     Model the variance of log-expression profiles for each gene,
#'     using a Poisson-based trend to estimate technical variance.
#'     Uses SQL aggregation for computing per-gene means and variances.
#'     \describe{
#'       \item{\code{size.factors}}{numeric vector of size factors for normalization}
#'       \item{\code{block}}{factor specifying blocking levels for each cell}
#'       \item{\code{design}}{design matrix for blocking (not supported with block)}
#'       \item{\code{subset.row}}{rows for which to model the variance}
#'       \item{\code{equiweight}}{whether blocks should be weighted equally}
#'       \item{\code{method}}{method for combining p-values across blocks}
#'     }
#'   }
#'   \item{\code{modelGeneCV2(x, size.factors = NULL, block = NULL,
#'     subset.row = NULL, subset.fit = NULL, ..., equiweight = TRUE, method = "fisher",
#'     BPPARAM = SerialParam())}:}{
#'     Model the squared coefficient of variation (CV2) of normalized expression
#'     for each gene, fitting a trend to account for the mean-CV2 relationship.
#'     Uses SQL aggregation for computing per-gene means and variances on
#'     size-factor normalized counts.
#'     \describe{
#'       \item{\code{size.factors}}{numeric vector of size factors for normalization}
#'       \item{\code{block}}{factor specifying blocking levels for each cell}
#'       \item{\code{subset.row}}{rows for which to model the CV2}
#'       \item{\code{subset.fit}}{rows to be used for trend fitting}
#'       \item{\code{equiweight}}{whether blocks should be weighted equally}
#'       \item{\code{method}}{method for combining p-values across blocks}
#'     }
#'   }
#' }
#'
#' @section Marker Gene Detection Methods:
#' The following marker gene detection methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{pairwiseTTests(x, groups, block = NULL, design = NULL,
#'     restrict = NULL, exclude = NULL, direction = c("any", "up", "down"),
#'     lfc = 0, std.lfc = FALSE, log.p = FALSE, gene.names = NULL,
#'     subset.row = NULL, BPPARAM = SerialParam())}:}{
#'     Perform pairwise Welch t-tests between groups of cells.
#'     Uses SQL aggregation for computing per-gene, per-group means and variances.
#'     \describe{
#'       \item{\code{groups}}{vector specifying group assignment for each cell}
#'       \item{\code{block}}{factor specifying blocking levels}
#'       \item{\code{restrict}}{character vector of groups to restrict comparisons to}
#'       \item{\code{exclude}}{character vector of groups to exclude from comparisons}
#'       \item{\code{direction}}{direction of log-fold changes to consider}
#'       \item{\code{lfc}}{log-fold change threshold}
#'       \item{\code{std.lfc}}{whether to standardize log-fold changes}
#'       \item{\code{log.p}}{whether to return log-transformed p-values}
#'       \item{\code{gene.names}}{character vector of gene names for output}
#'     }
#'   }
#'   \item{\code{pairwiseBinom(x, groups, block = NULL, restrict = NULL,
#'     exclude = NULL, direction = c("any", "up", "down"), threshold = 1e-8,
#'     lfc = 0, log.p = FALSE, gene.names = NULL, subset.row = NULL,
#'     BPPARAM = SerialParam())}:}{
#'     Perform pairwise binomial tests between groups of cells.
#'     Uses SQL aggregation for computing per-gene, per-group detection counts.
#'     \describe{
#'       \item{\code{groups}}{vector specifying group assignment for each cell}
#'       \item{\code{block}}{factor specifying blocking levels}
#'       \item{\code{restrict}}{character vector of groups to restrict comparisons to}
#'       \item{\code{exclude}}{character vector of groups to exclude from comparisons}
#'       \item{\code{direction}}{direction of log-fold changes to consider}
#'       \item{\code{threshold}}{expression threshold for detection}
#'       \item{\code{lfc}}{log-fold change threshold for proportions}
#'       \item{\code{log.p}}{whether to return log-transformed p-values}
#'       \item{\code{gene.names}}{character vector of gene names for output}
#'     }
#'   }
#'   \item{\code{findMarkers(x, groups, test.type = c("t", "wilcox", "binom"), ...,
#'     pval.type = c("any", "some", "all"), min.prop = NULL, log.p = FALSE,
#'     full.stats = FALSE, sorted = TRUE, row.data = NULL, add.summary = FALSE,
#'     BPPARAM = SerialParam())}:}{
#'     Find candidate marker genes for groups of cells.
#'     \describe{
#'       \item{\code{groups}}{vector specifying group assignment for each cell}
#'       \item{\code{test.type}}{type of pairwise test ("t" and "binom" are optimized)}
#'       \item{\code{pval.type}}{how to combine p-values across pairwise tests}
#'       \item{\code{sorted}}{whether to sort results by significance}
#'       \item{\code{add.summary}}{whether to add summary statistics}
#'     }
#'   }
#'   \item{\code{scoreMarkers(x, groups, block = NULL, pairings = NULL, lfc = 0,
#'     row.data = NULL, full.stats = FALSE, subset.row = NULL,
#'     BPPARAM = SerialParam(), true.auc = FALSE)}:}{
#'     Compute effect size summary statistics for potential marker genes.
#'     Uses SQL aggregation for computing per-gene, per-group means, variances,
#'     and detection proportions.
#'
#'     \strong{AUC computation options:} By default (\code{true.auc = FALSE}), AUC is
#'     computed using a normal approximation which is fast (~97\% correlation with scran).
#'     Set \code{true.auc = TRUE} to compute true rank-based AUC using SQL window functions
#'     (Wilcoxon rank-sum statistic). This is slower (~15-25x) but provides exact AUC
#'     values identical to scran (100\% correlation).
#'
#'     \describe{
#'       \item{\code{groups}}{vector specifying group assignment for each cell}
#'       \item{\code{block}}{factor specifying blocking levels}
#'       \item{\code{pairings}}{specification of pairwise comparisons to compute}
#'       \item{\code{lfc}}{log-fold change threshold for effect sizes}
#'       \item{\code{row.data}}{additional row data to include in output}
#'       \item{\code{full.stats}}{whether to return full pairwise statistics}
#'       \item{\code{true.auc}}{logical indicating whether to compute true rank-based
#'         AUC via SQL window functions (default FALSE uses normal approximation with
#'         ~97\% correlation; TRUE gives 100\% exact match but is ~15-25x slower)}
#'     }
#'   }
#'   \item{\code{summaryMarkerStats(x, groups, row.data = NULL, average = "mean",
#'     BPPARAM = SerialParam())}:}{
#'     Compute summary statistics for marker gene selection.
#'     Uses SQL aggregation for per-group means and detection proportions.
#'     \describe{
#'       \item{\code{groups}}{vector specifying group assignment for each cell}
#'       \item{\code{row.data}}{additional row data to include in output}
#'       \item{\code{average}}{type of average to compute ("mean" or "median")}
#'     }
#'   }
#' }
#'
#' @details
#' The DuckDBMatrix methods compute per-gene statistics using SQL
#' aggregation, which is efficient for large disk-backed matrices:
#' \itemize{
#'   \item For \code{pairwiseTTests}: per-group means and variances via SQL
#'   \item For \code{pairwiseBinom}: per-group detection counts via SQL
#' }
#' The trend fitting and p-value computations are then performed in R
#' using the standard scran functions.
#'
#' For \code{findMarkers}, both \code{test.type="t"} and \code{test.type="binom"}
#' use optimized DuckDBMatrix implementations. The Wilcoxon test
#' (\code{test.type="wilcox"}) falls back to the default method because it
#' requires pairwise cell comparisons that are not efficiently expressed in SQL.
#'
#' @section AUC Computation in scoreMarkers:
#' The \code{scoreMarkers} method for DuckDBMatrix provides two AUC computation modes:
#'
#' \strong{Normal approximation (default, \code{true.auc = FALSE}):}
#' \deqn{AUC \approx \Phi\left(\frac{\mu_1 - \mu_2}{\sqrt{\sigma_1^2/n_1 + \sigma_2^2/n_2}}\right)}
#' where \eqn{\Phi} is the standard normal CDF. This is fast but may be inaccurate for
#' sparse/bimodal data.
#'
#' \strong{True rank-based AUC (\code{true.auc = TRUE}):}
#' Computes exact AUC via the Wilcoxon rank-sum statistic using SQL window functions:
#' \deqn{AUC = \frac{R_1 - n_1(n_1+1)/2}{n_1 \cdot n_2}}
#' where \eqn{R_1} is the sum of ranks in group 1. This is ~15-25x slower than the
#' normal approximation but provides exact results identical to scran (correlation = 1.0).
#'
#' Benchmarks (5000 genes x 5000 cells, 12 pairwise comparisons):
#' \itemize{
#'   \item Default (normal approx): 0.68s, ~97\% correlation with scran
#'   \item \code{true.auc = TRUE}: 17.7s, 100\% exact match with scran
#' }
#'
#' The true AUC is preferred for:
#' \itemize{
#'   \item Genes with bimodal expression (expressed vs. not expressed)
#'   \item Sparse data with many zeros
#'   \item Applications where exact effect sizes are important
#' }
#'
#' Alternatively, use \code{findMarkers(x, groups, test.type = "wilcox")} for
#' rank-based p-values.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{DuckDBMatrix-class}} for the main class
#'   \item \code{\link{DuckDBMatrix-scuttle}} for the scuttle methods
#'   \item \code{\link[scran]{correlatePairs}} for the scran generic
#'   \item \code{\link[scran]{modelGeneVar}} for the scran generic
#'   \item \code{\link[scran]{modelGeneVarByPoisson}} for the scran generic
#'   \item \code{\link[scran]{modelGeneCV2}} for the scran generic
#'   \item \code{\link[scran]{pairwiseTTests}} for the scran generic
#'   \item \code{\link[scran]{pairwiseBinom}} for the scran generic
#'   \item \code{\link[scran]{findMarkers}} for the scran generic
#'   \item \code{\link[scran]{scoreMarkers}} for the scran generic
#'   \item \code{\link[scran]{summaryMarkerStats}} for the scran generic
#' }
#'
#' @aliases
#' correlatePairs,DuckDBMatrix-method
#' modelGeneVar,DuckDBMatrix-method
#' modelGeneVarByPoisson,DuckDBMatrix-method
#' modelGeneCV2,DuckDBMatrix-method
#' pairwiseTTests,DuckDBMatrix-method
#' pairwiseBinom,DuckDBMatrix-method
#' findMarkers,DuckDBMatrix-method
#' scoreMarkers,DuckDBMatrix-method
#' summaryMarkerStats,DuckDBMatrix-method
#'
#' @keywords utilities methods
#'
#' @name DuckDBMatrix-scran
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for computing per-group statistics
###

# Helper to register a data frame as a temporary DuckDB view (zero-copy)
# Returns the view name for use in joins
# OPTIMIZATION: Avoids copying data via left_join(..., copy = TRUE)
#' @importFrom duckdb duckdb_register duckdb_unregister
.register_temp_table <- function(db_conn, df, prefix = "temp") {
    # Generate unique table name to avoid collisions
    table_name <- paste0(prefix, "_", as.integer(Sys.time()), "_",
                         sample.int(10000, 1))
    duckdb::duckdb_register(db_conn, table_name, df, overwrite = TRUE)
    table_name
}

# Helper to unregister temporary tables
.unregister_temp_tables <- function(db_conn, table_names) {
    for (name in table_names) {
        try(duckdb::duckdb_unregister(db_conn, name), silent = TRUE)
    }
}

# Helper to extract table metadata from a DuckDBMatrix
# Consolidates repeated metadata extraction across helper functions
#' @importFrom DuckDBDataFrame keynames
.extract_table_metadata <- function(x) {
    table <- x@seed@table
    list(
        table = table,
        row_key = keynames(table)[1L],
        col_key = keynames(table)[2L],
        row_keycol = table@keycols[[1L]],
        col_keycol = table@keycols[[2L]],
        datacol = table@datacols[[1L]],
        fill = x@seed@fill
    )
}

# Helper to prepare block factor for blocked computations
# Consolidates repeated block factor preparation across helper functions
.prepare_block_factor <- function(block, ncells) {
    if (is.null(block)) {
        block <- factor(rep(1L, ncells))
        bnames <- NULL
    } else {
        if (length(block) != ncells) {
            stop("length of 'block' should be the same as 'ncol(x)'")
        }
        block <- as.factor(block)
        bnames <- levels(block)
    }
    nblocks <- nlevels(block)
    block_ncells <- as.integer(table(block))
    names(block_ncells) <- levels(block)
    
    list(
        block = block,
        bnames = bnames,
        nblocks = nblocks,
        block_ncells = block_ncells
    )
}

# Helper to prepare combined group × block factor
# Consolidates repeated group/block combination across helper functions
.prepare_group_block_factor <- function(groups, block, ncells) {
    groups <- as.factor(groups)
    group_levels <- levels(groups)
    ngroups <- nlevels(groups)
    
    if (is.null(block)) {
        combined <- groups
        combined_levels <- group_levels
        nblocks <- 1L
    } else {
        if (length(block) != ncells) {
            stop("length of 'block' does not equal 'ncol(x)'")
        }
        block <- as.factor(block)
        nblocks <- nlevels(block)
        combined <- factor(
            paste0(as.integer(groups), "_", as.integer(block)),
            levels = as.character(outer(seq_len(ngroups), seq_len(nblocks),
                                        paste, sep = "_"))
        )
        combined_levels <- levels(combined)
    }
    
    # Compute group sizes (excluding NAs)
    valid <- !is.na(groups)
    group_sizes <- as.integer(table(combined[valid]))
    names(group_sizes) <- combined_levels
    
    list(
        combined = combined,
        combined_levels = combined_levels,
        group_levels = group_levels,
        ngroups = ngroups,
        nblocks = nblocks,
        group_sizes = group_sizes,
        valid = valid
    )
}

# Compute per-gene means and variances, optionally by block
# Uses numerically stable two-pass algorithm:
#   Pass 1: Compute means
#   Pass 2: Compute sum of squared deviations from mean
# OPTIMIZATION: Uses duckdb_register for zero-copy joins instead of copy = TRUE
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter group_by left_join mutate summarize tbl
#' @importFrom stats setNames
.compute_blocked_stats_DuckDBMatrix <-
function(x, block = NULL, subset.row = NULL)
{
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    ddb_table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol
    fill <- meta$fill

    # Prepare block factor (REFACTORED)
    block_info <- .prepare_block_factor(block, ncells)
    block <- block_info$block
    bnames <- block_info$bnames
    nblocks <- block_info$nblocks
    block_ncells <- block_info$block_ncells

    # Create block assignment data frame
    block_df <- data.frame(
        col_idx = unname(col_keycol),
        block_id = as.character(block),
        stringsAsFactors = FALSE
    )
    names(block_df)[1L] <- col_key

    conn <- tblconn(ddb_table, select = FALSE)
    db_conn <- conn$src$con
    row_sym <- as.name(row_key)
    block_sym <- as.name("block_id")

    # Register block_df as temporary table (OPTIMIZATION: zero-copy)
    block_tbl_name <- .register_temp_table(db_conn, block_df, "block")
    temp_tables <- block_tbl_name
    on.exit(.unregister_temp_tables(db_conn, temp_tables), add = TRUE)

    # Pass 1: Compute means per gene per block
    mean_aggr <- list(
        stored_sum = call("sum", datacol, na.rm = TRUE),
        stored_n = call("n")
    )

    mean_result <- conn
    mean_result <- left_join(mean_result, tbl(db_conn, block_tbl_name), by = col_key)
    mean_result <- group_by(mean_result, !!row_sym, !!block_sym)
    mean_result <- summarize(mean_result, !!!mean_aggr, .groups = "drop")
    mean_result <- collect(mean_result)
    mean_result$stored_n <- as.integer(mean_result$stored_n)

    # Pivot to matrices and compute means
    row_keyvals <- unname(row_keycol)
    means <- matrix(0, nrow = ngenes, ncol = nblocks,
                    dimnames = list(rnames, bnames))
    vars <- matrix(0, nrow = ngenes, ncol = nblocks,
                   dimnames = list(rnames, bnames))

    # First compute all means (needed for pass 2)
    gene_means_by_block <- list()
    for (b in seq_len(nblocks)) {
        block_name <- levels(block)[b]
        block_n <- block_ncells[b]
        block_result <- mean_result[mean_result$block_id == block_name, , drop = FALSE]
        row_idx <- match(block_result[[row_key]], row_keyvals)

        stored_sum <- rep(fill * block_n, ngenes)
        stored_sum[row_idx] <- block_result$stored_sum +
            fill * (block_n - block_result$stored_n)
        means[, b] <- stored_sum / block_n
        gene_means_by_block[[block_name]] <- means[, b]
    }

    # Pass 2: Compute sum of squared deviations using two-pass algorithm
    # Create mean lookup table for SQL join
    mean_lookup <- data.frame(
        row_key = rep(row_keyvals, nblocks),
        block_id = rep(levels(block), each = ngenes),
        gene_mean = as.vector(means),
        stringsAsFactors = FALSE
    )
    names(mean_lookup)[1L] <- row_key

    # Register mean_lookup as temporary table (OPTIMIZATION: zero-copy)
    mean_tbl_name <- .register_temp_table(db_conn, mean_lookup, "means")
    temp_tables <- c(temp_tables, mean_tbl_name)

    # SQL: compute sum of (x - mean)^2 per gene per block
    dev_sq_sym <- as.name("dev_sq")
    mean_sym <- as.name("gene_mean")

    conn2 <- tblconn(ddb_table, select = FALSE)
    conn2 <- left_join(conn2, tbl(db_conn, block_tbl_name), by = col_key)
    conn2 <- left_join(conn2, tbl(db_conn, mean_tbl_name), by = c(row_key, "block_id"))

    # Compute (x - mean)^2
    dev_sq_expr <- list(dev_sq = call("*",
                                      call("-", datacol, mean_sym),
                                      call("-", datacol, mean_sym)))
    conn2 <- mutate(conn2, !!!dev_sq_expr)

    # Sum squared deviations per gene per block
    var_aggr <- list(
        sum_dev_sq = call("sum", dev_sq_sym, na.rm = TRUE),
        stored_n = call("n")
    )
    var_result <- group_by(conn2, !!row_sym, !!block_sym)
    var_result <- summarize(var_result, !!!var_aggr, .groups = "drop")
    var_result <- collect(var_result)
    var_result$stored_n <- as.integer(var_result$stored_n)

    # Compute variances with fill value contribution
    for (b in seq_len(nblocks)) {
        block_name <- levels(block)[b]
        block_n <- block_ncells[b]

        if (block_n <= 1L) {
            vars[, b] <- NA_real_
            next
        }

        block_var_result <- var_result[var_result$block_id == block_name, , drop = FALSE]
        row_idx <- match(block_var_result[[row_key]], row_keyvals)

        # Initialize with fill value contribution
        # For genes not in result, all values are fill, so variance contribution is 0
        # (fill - mean)^2 * n, but mean = fill when all values are fill
        gene_mean <- means[, b]
        n_fill <- block_n  # Start with all fill values
        sum_dev_sq <- (fill - gene_mean)^2 * n_fill

        # Update with actual data
        if (nrow(block_var_result) > 0L) {
            actual_n <- block_var_result$stored_n
            actual_sum_dev_sq <- block_var_result$sum_dev_sq
            n_fill_actual <- block_n - actual_n

            # Total sum of squared deviations =
            #   sum_dev_sq from actual values + (fill - mean)^2 * n_fill
            sum_dev_sq[row_idx] <- actual_sum_dev_sq +
                (fill - gene_mean[row_idx])^2 * n_fill_actual
        }

        # Sample variance with Bessel's correction
        vars[, b] <- sum_dev_sq / (block_n - 1L)
    }

    list(means = means, vars = vars, ncells = block_ncells)
}

# Compute per-gene means and variances by group (for t-tests)
# Uses numerically stable two-pass algorithm:
#   Pass 1: Compute means
#   Pass 2: Compute sum of squared deviations from mean
# OPTIMIZATION: Uses duckdb_register for zero-copy joins instead of copy = TRUE
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter group_by left_join mutate summarize tbl
#' @importFrom stats setNames
.compute_group_stats_DuckDBMatrix <-
function(x, groups, block = NULL, subset.row = NULL)
{
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    ddb_table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol
    fill <- meta$fill

    # Prepare group × block factor (REFACTORED)
    gb_info <- .prepare_group_block_factor(groups, block, ncells)
    combined <- gb_info$combined
    combined_levels <- gb_info$combined_levels
    group_levels <- gb_info$group_levels
    ngroups <- gb_info$ngroups
    nblocks <- gb_info$nblocks
    group_sizes <- gb_info$group_sizes
    valid <- gb_info$valid

    # Create group assignment data frame
    group_df <- data.frame(
        col_idx = unname(col_keycol),
        group_id = as.character(combined),
        stringsAsFactors = FALSE
    )
    names(group_df)[1L] <- col_key
    if (!all(valid)) {
        group_df <- group_df[valid, , drop = FALSE]
    }

    conn <- tblconn(ddb_table, select = FALSE)
    db_conn <- conn$src$con
    row_sym <- as.name(row_key)
    group_sym <- as.name("group_id")

    # Register group_df as temporary table (OPTIMIZATION: zero-copy)
    group_tbl_name <- .register_temp_table(db_conn, group_df, "group")
    temp_tables <- group_tbl_name
    on.exit(.unregister_temp_tables(db_conn, temp_tables), add = TRUE)

    # Pass 1: Compute means per gene per group
    mean_aggr <- list(
        stored_sum = call("sum", datacol, na.rm = TRUE),
        stored_n = call("n")
    )

    mean_result <- conn
    mean_result <- filter(mean_result, !!as.name(col_key) %in% local(group_df[[col_key]]))
    mean_result <- left_join(mean_result, tbl(db_conn, group_tbl_name), by = col_key)
    mean_result <- group_by(mean_result, !!row_sym, !!group_sym)
    mean_result <- summarize(mean_result, !!!mean_aggr, .groups = "drop")
    mean_result <- collect(mean_result)
    mean_result$stored_n <- as.integer(mean_result$stored_n)

    # Pivot to matrices and compute means
    row_keyvals <- unname(row_keycol)
    means <- matrix(NA_real_, nrow = ngenes, ncol = length(combined_levels),
                    dimnames = list(rnames, combined_levels))
    vars <- matrix(NA_real_, nrow = ngenes, ncol = length(combined_levels),
                   dimnames = list(rnames, combined_levels))

    # First compute all means (needed for pass 2)
    for (g in combined_levels) {
        group_n <- group_sizes[g]
        if (group_n == 0L) next

        group_mean_result <- mean_result[mean_result$group_id == g, , drop = FALSE]
        row_idx <- match(group_mean_result[[row_key]], row_keyvals)

        stored_sum <- rep(fill * group_n, ngenes)
        if (nrow(group_mean_result) > 0L) {
            stored_sum[row_idx] <- group_mean_result$stored_sum +
                fill * (group_n - group_mean_result$stored_n)
        }
        means[, g] <- stored_sum / group_n
    }

    # Pass 2: Compute sum of squared deviations using two-pass algorithm
    # Create mean lookup table for SQL join
    mean_lookup <- data.frame(
        row_key = rep(row_keyvals, length(combined_levels)),
        group_id = rep(combined_levels, each = ngenes),
        gene_mean = as.vector(means),
        stringsAsFactors = FALSE
    )
    names(mean_lookup)[1L] <- row_key

    # Register mean_lookup as temporary table (OPTIMIZATION: zero-copy)
    mean_tbl_name <- .register_temp_table(db_conn, mean_lookup, "means")
    temp_tables <- c(temp_tables, mean_tbl_name)

    # SQL: compute sum of (x - mean)^2 per gene per group
    dev_sq_sym <- as.name("dev_sq")
    mean_sym <- as.name("gene_mean")

    conn2 <- tblconn(ddb_table, select = FALSE)
    conn2 <- filter(conn2, !!as.name(col_key) %in% local(group_df[[col_key]]))
    conn2 <- left_join(conn2, tbl(db_conn, group_tbl_name), by = col_key)
    conn2 <- left_join(conn2, tbl(db_conn, mean_tbl_name), by = c(row_key, "group_id"))

    # Compute (x - mean)^2
    dev_sq_expr <- list(dev_sq = call("*",
                                      call("-", datacol, mean_sym),
                                      call("-", datacol, mean_sym)))
    conn2 <- mutate(conn2, !!!dev_sq_expr)

    # Sum squared deviations per gene per group
    var_aggr <- list(
        sum_dev_sq = call("sum", dev_sq_sym, na.rm = TRUE),
        stored_n = call("n")
    )
    var_result <- group_by(conn2, !!row_sym, !!group_sym)
    var_result <- summarize(var_result, !!!var_aggr, .groups = "drop")
    var_result <- collect(var_result)
    var_result$stored_n <- as.integer(var_result$stored_n)

    # Compute variances with fill value contribution
    for (g in combined_levels) {
        group_n <- group_sizes[g]
        if (group_n == 0L || group_n <= 1L) {
            vars[, g] <- NA_real_
            next
        }

        group_var_result <- var_result[var_result$group_id == g, , drop = FALSE]
        row_idx <- match(group_var_result[[row_key]], row_keyvals)

        # Initialize with fill value contribution
        gene_mean <- means[, g]
        n_fill <- group_n
        sum_dev_sq <- (fill - gene_mean)^2 * n_fill

        # Update with actual data
        if (nrow(group_var_result) > 0L) {
            actual_n <- group_var_result$stored_n
            actual_sum_dev_sq <- group_var_result$sum_dev_sq
            n_fill_actual <- group_n - actual_n

            sum_dev_sq[row_idx] <- actual_sum_dev_sq +
                (fill - gene_mean[row_idx])^2 * n_fill_actual
        }

        # Sample variance with Bessel's correction
        vars[, g] <- sum_dev_sq / (group_n - 1L)
    }

    if (is.null(block)) {
        list(means = means, vars = vars, ncells = group_sizes)
    } else {
        # Reshape into block structure
        list(means = means, vars = vars, ncells = group_sizes,
             ngroups = ngroups, nblocks = nblocks, group_levels = group_levels)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### modelGeneVar
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran modelGeneVar fitTrendVar combineBlocks
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
setMethod("modelGeneVar", "DuckDBMatrix",
function(x, block = NULL, design = NULL, subset.row = NULL, subset.fit = NULL,
         ..., equiweight = TRUE, method = "fisher", BPPARAM = SerialParam())
{
    if (!is.null(design)) {
        # Fall back to default method for design matrix
        return(callNextMethod())
    }

    # Compute mean and variance for subset.row
    x.stats <- .compute_blocked_stats_DuckDBMatrix(x, block = block,
                                                   subset.row = subset.row)

    # Compute mean and variance for subset.fit (if different)
    if (is.null(subset.fit)) {
        fit.stats <- x.stats
    } else {
        fit.stats <- .compute_blocked_stats_DuckDBMatrix(x, block = block,
                                                         subset.row = subset.fit)
    }

    # Decompose variance using scran's trend fitting
    collected <- .decompose_log_exprs_DuckDBMatrix(
        x.stats$means, x.stats$vars,
        fit.stats$means, fit.stats$vars,
        x.stats$ncells, ...
    )

    # Combine blocked statistics
    output <- combineBlocks(
        collected,
        method = method,
        equiweight = equiweight,
        weights = x.stats$ncells,
        valid = x.stats$ncells >= 2L,
        geometric = FALSE,
        ave.fields = c("mean", "total", "tech", "bio"),
        pval.field = "p.value"
    )

    # Set row names
    if (!is.null(subset.row)) {
        rownames(output) <- rownames(x)[subset.row]
    } else {
        rownames(output) <- rownames(x)
    }

    output
})

#' @importFrom scran fitTrendVar
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
.decompose_log_exprs_DuckDBMatrix <-
function(x.means, x.vars, fit.means, fit.vars, ncells, ...)
{
    dummy.trend.fit <- list(trend = function(x) rep(NA_real_, length(x)),
                            std.dev = NA_real_)

    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        fm <- fit.means[, i]
        fv <- fit.vars[, i]
        if (ncells[i] >= 2L) {
            fit <- fitTrendVar(fm, fv, ...)
        } else {
            fit <- dummy.trend.fit
        }

        xm <- unname(x.means[, i])
        xv <- unname(x.vars[, i])
        output <- DataFrame(mean = xm, total = xv, tech = fit$trend(xm))
        output$bio <- output$total - output$tech
        output$p.value <- pnorm(output$bio / output$tech, sd = fit$std.dev,
                                lower.tail = FALSE)
        output$FDR <- p.adjust(output$p.value, method = "BH")

        rownames(output) <- rownames(x.means)
        metadata(output) <- c(list(mean = fm, var = fv), fit)
        collected[[i]] <- output
    }
    names(collected) <- colnames(x.means)
    collected
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pairwiseTTests
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran pairwiseTTests combineMarkers
#' @importFrom S4Vectors DataFrame
#' @importFrom stats pt p.adjust
setMethod("pairwiseTTests", "DuckDBMatrix",
function(x, groups, block = NULL, design = NULL, restrict = NULL, exclude = NULL,
         direction = c("any", "up", "down"), lfc = 0, std.lfc = FALSE,
         log.p = FALSE, gene.names = NULL, subset.row = NULL,
         BPPARAM = SerialParam())
{
    if (!is.null(design)) {
        # Fall back to default method for design matrix
        return(callNextMethod())
    }

    direction <- match.arg(direction)

    # Setup groups
    groups <- .setup_groups_DuckDBMatrix(groups, x, restrict, exclude)
    group_levels <- levels(groups)
    ngroups <- nlevels(groups)

    # Compute group statistics
    stats <- .compute_group_stats_DuckDBMatrix(x, groups, block = block,
                                               subset.row = subset.row)

    # Setup gene names
    if (is.null(gene.names)) {
        if (!is.null(subset.row)) {
            gene.names <- rownames(x)[subset.row]
        } else {
            gene.names <- rownames(x)
        }
        if (is.null(gene.names)) {
            gene.names <- as.character(seq_len(nrow(stats$means)))
        }
    }

    if (is.null(block)) {
        # No blocking: direct pairwise comparisons
        .pairwise_t_no_block(stats$means, stats$vars, stats$ncells,
                             group_levels, direction, lfc, std.lfc,
                             log.p, gene.names)
    } else {
        # With blocking: combine across blocks
        .pairwise_t_with_block(stats, group_levels, direction, lfc, std.lfc,
                               log.p, gene.names)
    }
})

.setup_groups_DuckDBMatrix <- function(groups, x, restrict, exclude)
{
    ncells <- ncol(x)
    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }

    if (!is.null(restrict)) {
        groups[!groups %in% restrict] <- NA
        if (is.factor(groups)) groups <- droplevels(groups)
    }
    if (!is.null(exclude)) {
        groups[groups %in% exclude] <- NA
        if (is.factor(groups)) groups <- droplevels(groups)
    }

    is.empty <- groups == ""
    if (any(is.empty, na.rm = TRUE)) {
        warning("replacing empty 'groups' with NA")
        groups[is.empty] <- NA
    }

    groups <- as.factor(groups)
    if (nlevels(groups) < 2L) {
        stop("need at least two unique levels in 'groups'")
    }
    groups
}

#' @importFrom S4Vectors DataFrame
#' @importFrom stats pt
.pairwise_t_no_block <-
function(means, vars, ncells, group_levels, direction, lfc, std.lfc,
         log.p, gene.names)
{
    ngroups <- length(group_levels)
    ngenes <- nrow(means)

    collected.stats <- list()
    collected.pairs <- list()
    counter <- 1L

    for (h in seq_len(ngroups)) {
        host <- group_levels[h]
        host.n <- ncells[host]
        host.mean <- means[, host]
        host.s2 <- vars[, host]

        for (t in seq_len(h - 1L)) {
            target <- group_levels[t]
            target.n <- ncells[target]
            target.mean <- means[, target]
            target.s2 <- vars[, target]

            # Welch t-test
            t.stats <- .get_t_test_stats_DuckDBMatrix(host.s2, target.s2,
                                                      host.n, target.n)

            cur.lfc <- host.mean - target.mean
            p.out <- .run_t_test_DuckDBMatrix(cur.lfc, t.stats$err,
                                              t.stats$df, lfc)

            hvt.p <- .choose_pvalues_DuckDBMatrix(p.out$left, p.out$right,
                                                  direction)
            tvh.p <- .choose_pvalues_DuckDBMatrix(p.out$right, p.out$left,
                                                  direction)

            effect.hvt <- cur.lfc
            effect.tvh <- -cur.lfc
            if (std.lfc) {
                # Cohen's d
                pooled.s2 <- ((host.n - 1) * host.s2 + (target.n - 1) * target.s2) /
                    (target.n + host.n - 2)
                is.zero <- cur.lfc == 0
                effect.hvt <- cur.lfc / sqrt(pooled.s2)
                effect.hvt[is.zero] <- 0
                effect.tvh <- -effect.hvt
            }

            collected.stats[[counter]] <- list(
                .create_stats_df(effect.hvt, hvt.p, gene.names, log.p),
                .create_stats_df(effect.tvh, tvh.p, gene.names, log.p)
            )
            collected.pairs[[counter]] <- DataFrame(first = c(host, target),
                                                    second = c(target, host))
            counter <- counter + 1L
        }
    }

    output <- list(
        statistics = unlist(collected.stats, recursive = FALSE),
        pairs = do.call(rbind, collected.pairs)
    )
    .reorder_output(output, group_levels)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom metapod combineParallelPValues averageParallelStats
.pairwise_t_with_block <-
function(stats, group_levels, direction, lfc, std.lfc, log.p, gene.names)
{
    ngroups <- stats$ngroups
    nblocks <- stats$nblocks
    ngenes <- nrow(stats$means)

    collected.stats <- list()
    collected.pairs <- list()
    counter <- 1L

    for (h in seq_len(ngroups)) {
        host <- group_levels[h]

        for (t in seq_len(h - 1L)) {
            target <- group_levels[t]

            all.forward <- all.reverse <- all.left <- all.right <- vector("list", nblocks)
            all.weight <- numeric(nblocks)
            valid.test <- logical(nblocks)

            for (b in seq_len(nblocks)) {
                # Get stats for this block
                host.key <- paste0(h, "_", b)
                target.key <- paste0(t, "_", b)

                host.n <- stats$ncells[host.key]
                target.n <- stats$ncells[target.key]

                if (is.na(host.n) || is.na(target.n) ||
                    host.n < 2L || target.n < 2L) {
                    valid.test[b] <- FALSE
                    all.weight[b] <- 0
                    next
                }

                host.mean <- stats$means[, host.key]
                host.s2 <- stats$vars[, host.key]
                target.mean <- stats$means[, target.key]
                target.s2 <- stats$vars[, target.key]

                t.stats <- .get_t_test_stats_DuckDBMatrix(host.s2, target.s2,
                                                          host.n, target.n)

                cur.lfc <- host.mean - target.mean
                p.out <- .run_t_test_DuckDBMatrix(cur.lfc, t.stats$err,
                                                  t.stats$df, lfc)

                effect <- cur.lfc
                if (std.lfc) {
                    pooled.s2 <- ((host.n - 1) * host.s2 +
                                      (target.n - 1) * target.s2) /
                        (target.n + host.n - 2)
                    is.zero <- cur.lfc == 0
                    effect <- cur.lfc / sqrt(pooled.s2)
                    effect[is.zero] <- 0
                }

                all.forward[[b]] <- effect
                all.reverse[[b]] <- -effect
                all.left[[b]] <- p.out$left
                all.right[[b]] <- p.out$right
                all.weight[b] <- 1 / (1 / host.n + 1 / target.n)
                valid.test[b] <- all(!is.na(t.stats$df))
            }

            # Combine across blocks
            if (any(valid.test)) {
                w <- all.weight[valid.test]
                com.left <- combineParallelPValues(all.left[valid.test],
                                                   method = "stouffer",
                                                   weights = w, log.p = TRUE)$p.value
                com.right <- combineParallelPValues(all.right[valid.test],
                                                    method = "stouffer",
                                                    weights = w, log.p = TRUE)$p.value

                hvt.p <- .choose_pvalues_DuckDBMatrix(com.left, com.right, direction)
                tvh.p <- .choose_pvalues_DuckDBMatrix(com.right, com.left, direction)

                forward.effect <- averageParallelStats(all.forward[valid.test], w)
                reverse.effect <- averageParallelStats(all.reverse[valid.test], w)
            } else {
                hvt.p <- tvh.p <- forward.effect <- reverse.effect <-
                    rep(NA_real_, ngenes)
                warning(paste("no within-block comparison between", host,
                              "and", target))
            }

            collected.stats[[counter]] <- list(
                .create_stats_df(forward.effect, hvt.p, gene.names, log.p,
                                 effect.name = "logFC"),
                .create_stats_df(reverse.effect, tvh.p, gene.names, log.p,
                                 effect.name = "logFC")
            )
            collected.pairs[[counter]] <- DataFrame(first = c(host, target),
                                                    second = c(target, host))
            counter <- counter + 1L
        }
    }

    output <- list(
        statistics = unlist(collected.stats, recursive = FALSE),
        pairs = do.call(rbind, collected.pairs)
    )
    .reorder_output(output, group_levels)
}

.get_t_test_stats_DuckDBMatrix <- function(host.s2, target.s2, host.n, target.n)
{
    host.df <- max(0L, host.n - 1L)
    target.df <- max(0L, target.n - 1L)

    host.s2 <- pmax(host.s2, 1e-8)
    target.s2 <- pmax(target.s2, 1e-8)

    if (host.df > 0L && target.df > 0L) {
        host.err <- host.s2 / host.n
        target.err <- target.s2 / target.n
        cur.err <- host.err + target.err
        cur.df <- cur.err^2 / (host.err^2 / host.df + target.err^2 / target.df)
    } else {
        cur.err <- cur.df <- NA_real_
    }
    list(err = cur.err, df = cur.df)
}

#' @importFrom stats pt
.run_t_test_DuckDBMatrix <- function(cur.lfc, cur.err, cur.df, thresh.lfc = 0)
{
    thresh.lfc <- abs(thresh.lfc)
    if (thresh.lfc == 0) {
        cur.t <- cur.lfc / sqrt(cur.err)
        left <- pt(cur.t, df = cur.df, lower.tail = TRUE, log.p = TRUE)
        right <- pt(cur.t, df = cur.df, lower.tail = FALSE, log.p = TRUE)
    } else {
        lower.t <- (cur.lfc + thresh.lfc) / sqrt(cur.err)
        left <- pt(lower.t, df = cur.df, lower.tail = TRUE, log.p = TRUE)

        upper.t <- (cur.lfc - thresh.lfc) / sqrt(cur.err)
        right <- pt(upper.t, df = cur.df, lower.tail = FALSE, log.p = TRUE)
    }
    list(left = left, right = right)
}

.choose_pvalues_DuckDBMatrix <- function(left, right, direction)
{
    if (direction == "up") {
        right
    } else if (direction == "down") {
        left
    } else {
        log.p.out <- pmin(left, right) + log(2)
        pmin(0, log.p.out)
    }
}

#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
.create_stats_df <- function(effect, p, gene.names, log.p, effect.name = "logFC")
{
    # Remove names from vectors to match standard scran output
    p <- as.vector(p)
    effect <- unname(effect)
    effect.list <- list(effect)
    names(effect.list) <- effect.name

    if (log.p) {
        DataFrame(effect.list, log.p.value = p, log.FDR = .logBH_DuckDBMatrix(p),
                  check.names = FALSE, row.names = gene.names)
    } else {
        DataFrame(effect.list, p.value = exp(p),
                  FDR = exp(.logBH_DuckDBMatrix(p)),
                  check.names = FALSE, row.names = gene.names)
    }
}

.logBH_DuckDBMatrix <- function(log.p.val)
{
    o <- order(log.p.val)
    repval <- log.p.val[o] + log(length(o) / seq_along(o))
    repval <- rev(cummin(rev(repval)))
    repval[o] <- repval
    repval
}

.reorder_output <- function(output, levels)
{
    o <- order(
        match(output$pairs$first, levels),
        match(output$pairs$second, levels)
    )
    output$statistics <- output$statistics[o]
    output$pairs <- output$pairs[o, ]
    output
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pairwiseBinom
###

# Compute per-gene detection counts by group (for binomial tests)
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter group_by left_join summarize
#' @importFrom stats setNames
.compute_group_detection_stats_DuckDBMatrix <-
function(x, groups, block = NULL, subset.row = NULL, threshold = 1e-8)
{
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    ddb_table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol
    fill <- meta$fill

    # Prepare group × block factor (REFACTORED)
    gb_info <- .prepare_group_block_factor(groups, block, ncells)
    combined <- gb_info$combined
    combined_levels <- gb_info$combined_levels
    group_levels <- gb_info$group_levels
    ngroups <- gb_info$ngroups
    nblocks <- gb_info$nblocks
    group_sizes <- gb_info$group_sizes
    valid <- gb_info$valid

    # Create group assignment data frame
    group_df <- data.frame(
        col_idx = unname(col_keycol),
        group_id = as.character(combined),
        stringsAsFactors = FALSE
    )
    names(group_df)[1L] <- col_key
    if (!all(valid)) {
        group_df <- group_df[valid, , drop = FALSE]
    }

    conn <- tblconn(ddb_table, select = FALSE)

    # Compute count of cells above threshold per gene per group
    row_sym <- as.name(row_key)
    group_sym <- as.name("group_id")

    # Use countif-like aggregation: count cells where value > threshold
    aggr <- list(
        n_detected = call("sum", call("as.integer", call(">", datacol, threshold)),
                          na.rm = TRUE),
        stored_n = call("n")
    )

    result <- conn
    result <- filter(result, !!as.name(col_key) %in% local(group_df[[col_key]]))
    result <- left_join(result, group_df, by = col_key, copy = TRUE)
    result <- group_by(result, !!row_sym, !!group_sym)
    result <- summarize(result, !!!aggr, .groups = "drop")
    result <- collect(result)

    # Convert integer64 to integer
    result$stored_n <- as.integer(result$stored_n)
    result$n_detected <- as.integer(result$n_detected)

    # Pivot to matrices
    row_keyvals <- unname(row_keycol)

    # Create output matrix for detection counts
    ndetected <- matrix(0L, nrow = ngenes, ncol = length(combined_levels),
                        dimnames = list(rnames, combined_levels))

    # Handle fill value for genes not in result
    fill_detected <- as.integer(fill > threshold)

    for (g in combined_levels) {
        group_n <- group_sizes[g]
        if (group_n == 0L) next

        group_result <- result[result$group_id == g, , drop = FALSE]

        # Match rows
        row_idx <- match(group_result[[row_key]], row_keyvals)

        # Initialize with fill-based detection
        detected <- rep(fill_detected * group_n, ngenes)

        if (nrow(group_result) > 0L) {
            # Stored values detected + fill values detected for missing cells
            detected[row_idx] <- group_result$n_detected +
                fill_detected * (group_n - group_result$stored_n)
        }

        ndetected[, g] <- detected
    }

    if (is.null(block)) {
        list(ndetected = ndetected, ncells = group_sizes)
    } else {
        list(ndetected = ndetected, ncells = group_sizes, ngroups = ngroups,
             nblocks = nblocks, group_levels = group_levels)
    }
}

# Compute per-gene means, variances, AND detection counts by group in fewer SQL queries
# Optimization: Combines stats + detection into 2 SQL queries instead of 4
# Uses numerically stable two-pass algorithm for variance:
#   Pass 1: Compute sums AND detection counts together
#   Pass 2: Compute sum of squared deviations from mean
# OPTIMIZATION: Uses duckdb_register for zero-copy joins instead of copy = TRUE
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter group_by left_join mutate summarize tbl
#' @importFrom stats setNames
.compute_group_stats_with_detection_DuckDBMatrix <-
function(x, groups, block = NULL, subset.row = NULL, threshold = 1e-8)
{
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    ddb_table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol
    fill <- meta$fill

    # Prepare group × block factor (REFACTORED)
    gb_info <- .prepare_group_block_factor(groups, block, ncells)
    combined <- gb_info$combined
    combined_levels <- gb_info$combined_levels
    group_levels <- gb_info$group_levels
    ngroups <- gb_info$ngroups
    nblocks <- gb_info$nblocks
    group_sizes <- gb_info$group_sizes
    valid <- gb_info$valid

    # Create group assignment data frame
    group_df <- data.frame(
        col_idx = unname(col_keycol),
        group_id = as.character(combined),
        stringsAsFactors = FALSE
    )
    names(group_df)[1L] <- col_key
    if (!all(valid)) {
        group_df <- group_df[valid, , drop = FALSE]
    }

    conn <- tblconn(ddb_table, select = FALSE)
    db_conn <- conn$src$con
    row_sym <- as.name(row_key)
    group_sym <- as.name("group_id")

    # Register group_df as temporary table (OPTIMIZATION: zero-copy)
    group_tbl_name <- .register_temp_table(db_conn, group_df, "group")
    temp_tables <- group_tbl_name
    on.exit(.unregister_temp_tables(db_conn, temp_tables), add = TRUE)

    # Pass 1: Compute sums AND detection counts together (OPTIMIZED)
    # Previously this was 2 separate queries
    pass1_aggr <- list(
        stored_sum = call("sum", datacol, na.rm = TRUE),
        stored_n = call("n"),
        n_detected = call("sum", call("as.integer", call(">", datacol, threshold)),
                          na.rm = TRUE)
    )

    pass1_result <- conn
    pass1_result <- filter(pass1_result, !!as.name(col_key) %in% local(group_df[[col_key]]))
    pass1_result <- left_join(pass1_result, tbl(db_conn, group_tbl_name), by = col_key)
    pass1_result <- group_by(pass1_result, !!row_sym, !!group_sym)
    pass1_result <- summarize(pass1_result, !!!pass1_aggr, .groups = "drop")
    pass1_result <- collect(pass1_result)
    pass1_result$stored_n <- as.integer(pass1_result$stored_n)
    pass1_result$n_detected <- as.integer(pass1_result$n_detected)

    # Pivot to matrices and compute means
    row_keyvals <- unname(row_keycol)
    means <- matrix(NA_real_, nrow = ngenes, ncol = length(combined_levels),
                    dimnames = list(rnames, combined_levels))
    vars <- matrix(NA_real_, nrow = ngenes, ncol = length(combined_levels),
                   dimnames = list(rnames, combined_levels))
    ndetected <- matrix(0L, nrow = ngenes, ncol = length(combined_levels),
                        dimnames = list(rnames, combined_levels))

    # Handle fill value for detection
    fill_detected <- as.integer(fill > threshold)

    # Compute means and detection counts from pass 1 results
    for (g in combined_levels) {
        group_n <- group_sizes[g]
        if (group_n == 0L) next

        group_result <- pass1_result[pass1_result$group_id == g, , drop = FALSE]
        row_idx <- match(group_result[[row_key]], row_keyvals)

        # Initialize with fill values
        stored_sum <- rep(fill * group_n, ngenes)
        detected <- rep(fill_detected * group_n, ngenes)

        if (nrow(group_result) > 0L) {
            # Update sums with actual data
            stored_sum[row_idx] <- group_result$stored_sum +
                fill * (group_n - group_result$stored_n)
            # Update detection counts
            detected[row_idx] <- group_result$n_detected +
                fill_detected * (group_n - group_result$stored_n)
        }

        means[, g] <- stored_sum / group_n
        ndetected[, g] <- detected
    }

    # Pass 2: Compute sum of squared deviations (numerically stable)
    # Create mean lookup table for SQL join
    mean_lookup <- data.frame(
        row_key = rep(row_keyvals, length(combined_levels)),
        group_id = rep(combined_levels, each = ngenes),
        gene_mean = as.vector(means),
        stringsAsFactors = FALSE
    )
    names(mean_lookup)[1L] <- row_key

    # Register mean_lookup as temporary table (OPTIMIZATION: zero-copy)
    mean_tbl_name <- .register_temp_table(db_conn, mean_lookup, "means")
    temp_tables <- c(temp_tables, mean_tbl_name)

    # SQL: compute sum of (x - mean)^2 per gene per group
    dev_sq_sym <- as.name("dev_sq")
    mean_sym <- as.name("gene_mean")

    conn2 <- tblconn(ddb_table, select = FALSE)
    conn2 <- filter(conn2, !!as.name(col_key) %in% local(group_df[[col_key]]))
    conn2 <- left_join(conn2, tbl(db_conn, group_tbl_name), by = col_key)
    conn2 <- left_join(conn2, tbl(db_conn, mean_tbl_name), by = c(row_key, "group_id"))

    # Compute (x - mean)^2
    dev_sq_expr <- list(dev_sq = call("*",
                                      call("-", datacol, mean_sym),
                                      call("-", datacol, mean_sym)))
    conn2 <- mutate(conn2, !!!dev_sq_expr)

    # Sum squared deviations per gene per group
    var_aggr <- list(
        sum_dev_sq = call("sum", dev_sq_sym, na.rm = TRUE),
        stored_n = call("n")
    )
    var_result <- group_by(conn2, !!row_sym, !!group_sym)
    var_result <- summarize(var_result, !!!var_aggr, .groups = "drop")
    var_result <- collect(var_result)
    var_result$stored_n <- as.integer(var_result$stored_n)

    # Compute variances with fill value contribution
    for (g in combined_levels) {
        group_n <- group_sizes[g]
        if (group_n == 0L || group_n <= 1L) {
            vars[, g] <- NA_real_
            next
        }

        group_var_result <- var_result[var_result$group_id == g, , drop = FALSE]
        row_idx <- match(group_var_result[[row_key]], row_keyvals)

        # Initialize with fill value contribution
        gene_mean <- means[, g]
        n_fill <- group_n
        sum_dev_sq <- (fill - gene_mean)^2 * n_fill

        # Update with actual data
        if (nrow(group_var_result) > 0L) {
            actual_n <- group_var_result$stored_n
            actual_sum_dev_sq <- group_var_result$sum_dev_sq
            n_fill_actual <- group_n - actual_n

            sum_dev_sq[row_idx] <- actual_sum_dev_sq +
                (fill - gene_mean[row_idx])^2 * n_fill_actual
        }

        # Sample variance with Bessel's correction
        vars[, g] <- sum_dev_sq / (group_n - 1L)
    }

    if (is.null(block)) {
        list(means = means, vars = vars, ndetected = ndetected, ncells = group_sizes)
    } else {
        list(means = means, vars = vars, ndetected = ndetected, ncells = group_sizes,
             ngroups = ngroups, nblocks = nblocks, group_levels = group_levels)
    }
}

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran pairwiseBinom
#' @importFrom S4Vectors DataFrame
#' @importFrom stats pbinom
setMethod("pairwiseBinom", "DuckDBMatrix",
function(x, groups, block = NULL, restrict = NULL, exclude = NULL,
         direction = c("any", "up", "down"), threshold = 1e-8, lfc = 0,
         log.p = FALSE, gene.names = NULL, subset.row = NULL,
         BPPARAM = SerialParam())
{
    direction <- match.arg(direction)

    # Setup groups
    groups <- .setup_groups_DuckDBMatrix(groups, x, restrict, exclude)
    group_levels <- levels(groups)
    ngroups <- nlevels(groups)

    # Compute detection statistics
    stats <- .compute_group_detection_stats_DuckDBMatrix(x, groups, block = block,
                                                         subset.row = subset.row,
                                                         threshold = threshold)

    # Setup gene names
    if (is.null(gene.names)) {
        if (!is.null(subset.row)) {
            gene.names <- rownames(x)[subset.row]
        } else {
            gene.names <- rownames(x)
        }
        if (is.null(gene.names)) {
            gene.names <- as.character(seq_len(nrow(stats$ndetected)))
        }
    }

    if (is.null(block)) {
        .pairwise_binom_no_block(stats$ndetected, stats$ncells,
                                 group_levels, direction, lfc, log.p, gene.names)
    } else {
        .pairwise_binom_with_block(stats, group_levels, direction, lfc,
                                   log.p, gene.names)
    }
})

#' @importFrom S4Vectors DataFrame
#' @importFrom stats pbinom
.pairwise_binom_no_block <-
function(ndetected, ncells, group_levels, direction, lfc, log.p, gene.names)
{
    ngroups <- length(group_levels)
    ngenes <- nrow(ndetected)

    collected.stats <- list()
    collected.pairs <- list()
    counter <- 1L

    for (h in seq_len(ngroups)) {
        host <- group_levels[h]
        host.n <- ncells[host]
        host.nzero <- ndetected[, host]

        for (t in seq_len(h - 1L)) {
            target <- group_levels[t]
            target.n <- ncells[target]
            target.nzero <- ndetected[, target]

            # Binomial test
            binom.out <- .run_binom_test_DuckDBMatrix(host.nzero, host.n,
                                                      target.nzero, target.n,
                                                      direction, lfc)

            hvt.p <- .choose_pvalues_DuckDBMatrix(binom.out$left, binom.out$right,
                                                  direction)
            tvh.p <- .choose_pvalues_DuckDBMatrix(binom.out$right, binom.out$left,
                                                  direction)

            effect.hvt <- binom.out$effect
            effect.tvh <- -binom.out$effect

            collected.stats[[counter]] <- list(
                .create_stats_df(effect.hvt, hvt.p, gene.names, log.p,
                                 effect.name = "logFC"),
                .create_stats_df(effect.tvh, tvh.p, gene.names, log.p,
                                 effect.name = "logFC")
            )
            collected.pairs[[counter]] <- DataFrame(first = c(host, target),
                                                    second = c(target, host))
            counter <- counter + 1L
        }
    }

    output <- list(
        statistics = unlist(collected.stats, recursive = FALSE),
        pairs = do.call(rbind, collected.pairs)
    )
    .reorder_output(output, group_levels)
}

#' @importFrom S4Vectors DataFrame
#' @importFrom metapod combineParallelPValues averageParallelStats
.pairwise_binom_with_block <-
function(stats, group_levels, direction, lfc, log.p, gene.names)
{
    ngroups <- stats$ngroups
    nblocks <- stats$nblocks
    ngenes <- nrow(stats$ndetected)

    collected.stats <- list()
    collected.pairs <- list()
    counter <- 1L

    for (h in seq_len(ngroups)) {
        host <- group_levels[h]

        for (t in seq_len(h - 1L)) {
            target <- group_levels[t]

            # Collect per-block statistics
            all.left <- all.right <- all.effect <- all.weight <- vector("list", nblocks)
            all.valid <- logical(nblocks)

            for (b in seq_len(nblocks)) {
                combined_host <- paste0(h, "_", b)
                combined_target <- paste0(t, "_", b)

                host.n <- stats$ncells[combined_host]
                target.n <- stats$ncells[combined_target]

                if (is.na(host.n) || is.na(target.n)) {
                    host.n <- target.n <- 0L
                }

                all.valid[b] <- host.n > 0L && target.n > 0L
                all.weight[[b]] <- as.double(host.n) + as.double(target.n)

                if (all.valid[b]) {
                    host.nzero <- stats$ndetected[, combined_host]
                    target.nzero <- stats$ndetected[, combined_target]

                    binom.out <- .run_binom_test_DuckDBMatrix(host.nzero, host.n,
                                                              target.nzero, target.n,
                                                              direction, lfc)

                    all.left[[b]] <- binom.out$left
                    all.right[[b]] <- binom.out$right
                    all.effect[[b]] <- binom.out$effect
                } else {
                    all.left[[b]] <- all.right[[b]] <- all.effect[[b]] <- rep(NA_real_, ngenes)
                }
            }

            # Combine across blocks using Stouffer's method
            weights <- unlist(all.weight)
            valid_blocks <- which(all.valid)

            if (length(valid_blocks) == 0L) {
                hvt.p <- tvh.p <- rep(NA_real_, ngenes)
                effect.hvt <- effect.tvh <- rep(NA_real_, ngenes)
            } else {
                if (direction == "any") {
                    left.combined <- combineParallelPValues(all.left[valid_blocks],
                                                            method = "stouffer",
                                                            weights = weights[valid_blocks],
                                                            log.p = TRUE)$p.value
                    right.combined <- combineParallelPValues(all.right[valid_blocks],
                                                             method = "stouffer",
                                                             weights = weights[valid_blocks],
                                                             log.p = TRUE)$p.value
                    hvt.p <- .choose_pvalues_DuckDBMatrix(left.combined, right.combined,
                                                          direction)
                    tvh.p <- .choose_pvalues_DuckDBMatrix(right.combined, left.combined,
                                                          direction)
                } else if (direction == "up") {
                    right.combined <- combineParallelPValues(all.right[valid_blocks],
                                                             method = "stouffer",
                                                             weights = weights[valid_blocks],
                                                             log.p = TRUE)$p.value
                    hvt.p <- right.combined
                    tvh.p <- combineParallelPValues(all.left[valid_blocks],
                                                    method = "stouffer",
                                                    weights = weights[valid_blocks],
                                                    log.p = TRUE)$p.value
                } else {
                    left.combined <- combineParallelPValues(all.left[valid_blocks],
                                                            method = "stouffer",
                                                            weights = weights[valid_blocks],
                                                            log.p = TRUE)$p.value
                    hvt.p <- left.combined
                    tvh.p <- combineParallelPValues(all.right[valid_blocks],
                                                    method = "stouffer",
                                                    weights = weights[valid_blocks],
                                                    log.p = TRUE)$p.value
                }

                # Average effect sizes
                effect.hvt <- averageParallelStats(all.effect[valid_blocks],
                                                   weights = weights[valid_blocks])
                effect.tvh <- -effect.hvt
            }

            collected.stats[[counter]] <- list(
                .create_stats_df(effect.hvt, hvt.p, gene.names, log.p,
                                 effect.name = "logFC"),
                .create_stats_df(effect.tvh, tvh.p, gene.names, log.p,
                                 effect.name = "logFC")
            )
            collected.pairs[[counter]] <- DataFrame(first = c(host, target),
                                                    second = c(target, host))
            counter <- counter + 1L
        }
    }

    output <- list(
        statistics = unlist(collected.stats, recursive = FALSE),
        pairs = do.call(rbind, collected.pairs)
    )
    .reorder_output(output, group_levels)
}

#' @importFrom stats pbinom
.run_binom_test_DuckDBMatrix <-
function(host.nzero, host.n, target.nzero, target.n, direction, lfc)
{
    size <- host.nzero + target.nzero

    # Log-fold change in proportions, mimic edgeR::cpm()
    mean.lib <- mean(c(host.n, target.n))
    pseudo.host <- 1 * host.n / mean.lib
    pseudo.target <- 1 * target.n / mean.lib
    effect <- log2((host.nzero + pseudo.host) / (host.n + 2 * pseudo.host)) -
              log2((target.nzero + pseudo.target) / (target.n + 2 * pseudo.target))

    if (lfc == 0) {
        # Standard binomial test
        p <- host.n / (host.n + target.n)
        left <- pbinom(host.nzero, size, p, log.p = TRUE)
        right <- pbinom(host.nzero - 1, size, p, lower.tail = FALSE, log.p = TRUE)
    } else {
        # LFC threshold test
        fold <- 2^lfc
        p.left <- host.n / fold / (target.n + host.n / fold)
        p.right <- host.n * fold / (target.n + host.n * fold)
        left <- pbinom(host.nzero, size, p.left, log.p = TRUE)
        right <- pbinom(host.nzero - 1, size, p.right, lower.tail = FALSE, log.p = TRUE)
    }

    list(left = left, right = right, effect = effect)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findMarkers
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scran findMarkers combineMarkers summaryMarkerStats
#' @importFrom scuttle .bpNotSharedOrUp
setMethod("findMarkers", "DuckDBMatrix",
function(x, groups, test.type = c("t", "wilcox", "binom"), ...,
         pval.type = c("any", "some", "all"), min.prop = NULL, log.p = FALSE,
         full.stats = FALSE, sorted = TRUE, row.data = NULL,
         add.summary = FALSE, BPPARAM = SerialParam())
{
    test.type <- match.arg(test.type)

    # Only t-test and binom are optimized; wilcox falls back
    if (test.type == "wilcox") {
        return(callNextMethod())
    }

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Use optimized pairwise test
    if (test.type == "t") {
        fit <- pairwiseTTests(x, groups, ..., log.p = TRUE, BPPARAM = BPPARAM)
    } else {
        # test.type == "binom"
        fit <- pairwiseBinom(x, groups, ..., log.p = TRUE, BPPARAM = BPPARAM)
    }

    # Combine markers using scran's combineMarkers
    output <- combineMarkers(fit$statistics, fit$pairs, pval.type = pval.type,
                             min.prop = min.prop, log.p.in = TRUE,
                             log.p.out = log.p, full.stats = full.stats,
                             pval.field = "log.p.value", effect.field = "logFC",
                             sorted = sorted, BPPARAM = BPPARAM)

    # Add summary statistics if requested
    if (add.summary) {
        row.data <- summaryMarkerStats(as.matrix(x), groups, row.data = row.data,
                                       BPPARAM = BPPARAM)
    }

    # Add row data if provided
    if (!is.null(row.data)) {
        output <- .add_row_data_DuckDBMatrix(output, row.data,
                                              match.names = sorted)
    }

    output
})

#' @importClassesFrom S4Vectors DataFrame
#' @importFrom BiocGenerics cbind
.add_row_data_DuckDBMatrix <- function(output, row.data, match.names)
{
    if (is.null(row.data)) {
        return(output)
    }

    for (i in names(output)) {
        current <- output[[i]]

        if (is.data.frame(row.data) || is(row.data, "DataFrame")) {
            rd <- row.data
        } else {
            if (!i %in% names(row.data)) {
                stop("list-like 'row.data' should be named with the levels of 'groups'")
            }
            rd <- row.data[[i]]
        }

        rn <- rownames(current)
        if (match.names) {
            if (is.null(rn) || !identical(sort(rn), sort(rownames(rd)))) {
                stop("inconsistent or NULL row names for 'row.data' and result tables")
            }
            rd <- rd[rn, , drop = FALSE]
        } else if (!identical(rn, rownames(rd))) {
            stop("inconsistent row names for 'row.data' and result tables")
        }

        output[[i]] <- cbind(rd, current)
    }

    output
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### modelGeneVarByPoisson
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran modelGeneVarByPoisson fitTrendPoisson combineBlocks
#' @importFrom scuttle librarySizeFactors normalizeCounts
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
setMethod("modelGeneVarByPoisson", "DuckDBMatrix",
function(x, size.factors = NULL, block = NULL, design = NULL, subset.row = NULL,
         npts = 1000, dispersion = 0, pseudo.count = 1, ...,
         equiweight = TRUE, method = "fisher", BPPARAM = SerialParam())
{
    if (!is.null(design)) {
        # Fall back to default method for design matrix
        return(callNextMethod())
    }

    # Compute size factors if not provided
    if (is.null(size.factors)) {
        size.factors <- colSums(x)
        size.factors <- size.factors / mean(size.factors)
    }

    # Log-normalize the data before computing stats
    # modelGeneVarByPoisson expects log-normalized data for mean/var calculation
    log_x <- normalizeCounts(x, size.factors = size.factors,
                             log = TRUE, pseudo.count = pseudo.count)

    # Compute mean and variance on log-normalized data
    x.stats <- .compute_blocked_stats_DuckDBMatrix(log_x, block = block,
                                                   subset.row = subset.row)

    # Generate Poisson-based trend
    sim.out <- .generate_poisson_values_DuckDBMatrix(
        x.stats$means, size.factors, block = block, npts = npts,
        dispersion = dispersion, pseudo.count = pseudo.count
    )

    # Handle case where Poisson simulation failed
    if (length(sim.out$means) == 0L) {
        # Fall back to using actual data for trend
        sim.out <- list(
            means = as.vector(x.stats$means),
            vars = as.vector(x.stats$vars)
        )
    }

    # Decompose variance using Poisson trend
    collected <- .decompose_log_exprs_poisson_DuckDBMatrix(
        x.stats$means, x.stats$vars, sim.out$means, sim.out$vars,
        x.stats$ncells, ...
    )

    # Combine blocked statistics
    output <- combineBlocks(
        collected,
        method = method,
        equiweight = equiweight,
        weights = x.stats$ncells,
        valid = x.stats$ncells >= 2L,
        geometric = FALSE,
        ave.fields = c("mean", "total", "tech", "bio"),
        pval.field = "p.value"
    )

    # Set row names
    if (!is.null(subset.row)) {
        rownames(output) <- rownames(x)[subset.row]
    } else {
        rownames(output) <- rownames(x)
    }

    output
})

#' @importFrom BiocGenerics var
#' @importFrom stats rpois rnbinom
.generate_poisson_values_DuckDBMatrix <-
function(x.means, size.factors, block, npts, dispersion, pseudo.count)
{
    # Get range of means for simulation
    all_means <- as.vector(x.means)
    valid_means <- all_means[is.finite(all_means) & all_means > 0]
    if (length(valid_means) == 0L) {
        return(list(means = numeric(0), vars = numeric(0)))
    }

    # Compute limits in original count scale
    mean_range <- range(valid_means)
    xlim <- 2^mean_range - pseudo.count
    xlim[1L] <- max(xlim[1L], 1e-8)
    xlim[2L] <- min(xlim[2L], 1e6)  # Cap maximum to avoid Inf

    # Ensure finite limits
    if (!is.finite(xlim[1L]) || !is.finite(xlim[2L])) {
        xlim <- c(1e-8, 1e4)  # Use reasonable defaults
    }

    # Generate simulation points
    sim_means <- exp(seq(log(xlim[1L]), log(xlim[2L]), length.out = npts))

    # Simulate Poisson counts and compute log-normalized variance
    ncells <- length(size.factors)
    sim_vars <- numeric(npts)

    for (i in seq_len(npts)) {
        lambda <- sim_means[i] * size.factors
        if (dispersion > 0) {
            # Negative binomial
            counts <- rnbinom(ncells, mu = lambda, size = 1 / dispersion)
        } else {
            # Poisson
            counts <- rpois(ncells, lambda = lambda)
        }
        # Log-normalize
        log_vals <- log2(counts / size.factors + pseudo.count)
        sim_vars[i] <- var(log_vals)
    }

    # Convert back to log2 scale for means
    sim_log_means <- log2(sim_means + pseudo.count)

    list(means = sim_log_means, vars = sim_vars)
}

#' @importFrom scran fitTrendVar
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
.decompose_log_exprs_poisson_DuckDBMatrix <-
function(x.means, x.vars, fit.means, fit.vars, ncells, ...)
{
    # Fit trend to simulated Poisson data
    fit <- fitTrendVar(fit.means, fit.vars, ...)

    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        xm <- unname(x.means[, i])
        xv <- unname(x.vars[, i])

        output <- DataFrame(mean = xm, total = xv, tech = fit$trend(xm))
        output$bio <- output$total - output$tech
        output$p.value <- pnorm(output$bio / output$tech, sd = fit$std.dev,
                                lower.tail = FALSE)
        output$FDR <- p.adjust(output$p.value, method = "BH")

        rownames(output) <- rownames(x.means)
        metadata(output) <- c(list(mean = fit.means, var = fit.vars), fit)
        collected[[i]] <- output
    }
    names(collected) <- colnames(x.means)
    collected
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### modelGeneCV2
###

# Compute per-gene means and variances on size-factor normalized data
# Uses two-pass algorithm with size factor normalization incorporated in SQL
# OPTIMIZATION: Uses duckdb_register for zero-copy joins instead of copy = TRUE
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom dplyr collect filter group_by left_join mutate summarize tbl
#' @importFrom stats setNames
.compute_blocked_stats_normalized_DuckDBMatrix <-
function(x, size.factors, block = NULL, subset.row = NULL)
{
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ncells <- ncol(x)
    ngenes <- nrow(x)
    rnames <- rownames(x)

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    ddb_table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol
    fill <- meta$fill

    # Prepare block factor (REFACTORED)
    block_info <- .prepare_block_factor(block, ncells)
    block <- block_info$block
    bnames <- block_info$bnames
    nblocks <- block_info$nblocks
    block_ncells <- block_info$block_ncells

    # Create block and size factor assignment data frame
    sf_df <- data.frame(
        col_idx = unname(col_keycol),
        block_id = as.character(block),
        sf = size.factors,
        stringsAsFactors = FALSE
    )
    names(sf_df)[1L] <- col_key

    conn <- tblconn(ddb_table, select = FALSE)
    db_conn <- conn$src$con
    row_sym <- as.name(row_key)
    block_sym <- as.name("block_id")
    sf_sym <- as.name("sf")

    # Register sf_df as temporary table (OPTIMIZATION: zero-copy)
    sf_tbl_name <- .register_temp_table(db_conn, sf_df, "sf")
    temp_tables <- sf_tbl_name
    on.exit(.unregister_temp_tables(db_conn, temp_tables), add = TRUE)

    # Compute normalized value: value / sf
    norm_val_expr <- list(norm_val = call("/", datacol, sf_sym))

    # Pass 1: Compute means of normalized values per gene per block
    mean_aggr <- list(
        stored_sum = call("sum", as.name("norm_val"), na.rm = TRUE),
        stored_n = call("n")
    )

    mean_result <- conn
    mean_result <- left_join(mean_result, tbl(db_conn, sf_tbl_name), by = col_key)
    mean_result <- mutate(mean_result, !!!norm_val_expr)
    mean_result <- group_by(mean_result, !!row_sym, !!block_sym)
    mean_result <- summarize(mean_result, !!!mean_aggr, .groups = "drop")
    mean_result <- collect(mean_result)
    mean_result$stored_n <- as.integer(mean_result$stored_n)

    # Pivot to matrices and compute means
    row_keyvals <- unname(row_keycol)
    means <- matrix(0, nrow = ngenes, ncol = nblocks,
                    dimnames = list(rnames, bnames))
    vars <- matrix(0, nrow = ngenes, ncol = nblocks,
                   dimnames = list(rnames, bnames))

    # Normalized fill value (fill / mean(sf) approximately, but we use 0 for sparse)
    norm_fill <- fill / mean(size.factors)

    # First compute all means (needed for pass 2)
    for (b in seq_len(nblocks)) {
        block_name <- levels(block)[b]
        block_n <- block_ncells[b]
        block_result <- mean_result[mean_result$block_id == block_name, , drop = FALSE]
        row_idx <- match(block_result[[row_key]], row_keyvals)

        stored_sum <- rep(norm_fill * block_n, ngenes)
        if (nrow(block_result) > 0L) {
            stored_sum[row_idx] <- block_result$stored_sum +
                norm_fill * (block_n - block_result$stored_n)
        }
        means[, b] <- stored_sum / block_n
    }

    # Pass 2: Compute sum of squared deviations using two-pass algorithm
    # Create mean lookup table for SQL join
    mean_lookup <- data.frame(
        row_key = rep(row_keyvals, nblocks),
        block_id = rep(levels(block), each = ngenes),
        gene_mean = as.vector(means),
        stringsAsFactors = FALSE
    )
    names(mean_lookup)[1L] <- row_key

    # Register mean_lookup as temporary table (OPTIMIZATION: zero-copy)
    mean_tbl_name <- .register_temp_table(db_conn, mean_lookup, "means")
    temp_tables <- c(temp_tables, mean_tbl_name)

    # SQL: compute sum of (norm_val - mean)^2 per gene per block
    dev_sq_sym <- as.name("dev_sq")
    mean_sym <- as.name("gene_mean")
    norm_val_sym <- as.name("norm_val")

    conn2 <- tblconn(ddb_table, select = FALSE)
    conn2 <- left_join(conn2, tbl(db_conn, sf_tbl_name), by = col_key)
    conn2 <- mutate(conn2, !!!norm_val_expr)
    conn2 <- left_join(conn2, tbl(db_conn, mean_tbl_name), by = c(row_key, "block_id"))

    # Compute (norm_val - mean)^2
    dev_sq_expr <- list(dev_sq = call("*",
                                      call("-", norm_val_sym, mean_sym),
                                      call("-", norm_val_sym, mean_sym)))
    conn2 <- mutate(conn2, !!!dev_sq_expr)

    # Sum squared deviations per gene per block
    var_aggr <- list(
        sum_dev_sq = call("sum", dev_sq_sym, na.rm = TRUE),
        stored_n = call("n")
    )
    var_result <- group_by(conn2, !!row_sym, !!block_sym)
    var_result <- summarize(var_result, !!!var_aggr, .groups = "drop")
    var_result <- collect(var_result)
    var_result$stored_n <- as.integer(var_result$stored_n)

    # Compute variances with fill value contribution
    for (b in seq_len(nblocks)) {
        block_name <- levels(block)[b]
        block_n <- block_ncells[b]

        if (block_n <= 1L) {
            vars[, b] <- NA_real_
            next
        }

        block_var_result <- var_result[var_result$block_id == block_name, , drop = FALSE]
        row_idx <- match(block_var_result[[row_key]], row_keyvals)

        # Initialize with fill value contribution
        gene_mean <- means[, b]
        n_fill <- block_n
        sum_dev_sq <- (norm_fill - gene_mean)^2 * n_fill

        # Update with actual data
        if (nrow(block_var_result) > 0L) {
            actual_n <- block_var_result$stored_n
            actual_sum_dev_sq <- block_var_result$sum_dev_sq
            n_fill_actual <- block_n - actual_n

            sum_dev_sq[row_idx] <- actual_sum_dev_sq +
                (norm_fill - gene_mean[row_idx])^2 * n_fill_actual
        }

        # Sample variance with Bessel's correction
        vars[, b] <- sum_dev_sq / (block_n - 1L)
    }

    list(means = means, vars = vars, ncells = block_ncells)
}

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran modelGeneCV2 fitTrendCV2 combineBlocks
#' @importFrom scuttle librarySizeFactors
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
setMethod("modelGeneCV2", "DuckDBMatrix",
function(x, size.factors = NULL, block = NULL, subset.row = NULL,
         subset.fit = NULL, ..., equiweight = TRUE, method = "fisher",
         BPPARAM = SerialParam())
{
    # Compute size factors if not provided
    if (is.null(size.factors)) {
        if (!is.null(subset.row)) {
            size.factors <- colSums(x[subset.row, , drop = FALSE])
        } else {
            size.factors <- colSums(x)
        }
    }
    size.factors <- size.factors / mean(size.factors)

    # Compute mean and variance on size-factor normalized data via SQL
    x.stats <- .compute_blocked_stats_normalized_DuckDBMatrix(
        x, size.factors, block = block, subset.row = subset.row
    )

    # Compute stats for subset.fit (if different)
    if (is.null(subset.fit)) {
        fit.stats <- x.stats
    } else {
        fit.stats <- .compute_blocked_stats_normalized_DuckDBMatrix(
            x, size.factors, block = block, subset.row = subset.fit
        )
    }

    # Decompose CV2 using scran's trend fitting
    collected <- .decompose_cv2_DuckDBMatrix(
        x.stats$means, x.stats$vars,
        fit.stats$means, fit.stats$vars,
        x.stats$ncells, ...
    )

    # Combine blocked statistics
    output <- combineBlocks(
        collected,
        method = method,
        equiweight = equiweight,
        weights = x.stats$ncells,
        valid = x.stats$ncells >= 2L,
        geometric = TRUE,  # CV2 uses geometric mean for combining
        ave.fields = c("mean", "total", "trend", "ratio"),
        pval.field = "p.value"
    )

    # Set row names
    if (!is.null(subset.row)) {
        rownames(output) <- rownames(x)[subset.row]
    } else {
        rownames(output) <- rownames(x)
    }

    output
})

#' @importFrom scran fitTrendCV2
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom stats pnorm p.adjust
.decompose_cv2_DuckDBMatrix <-
function(x.means, x.vars, fit.means, fit.vars, ncells, ...)
{
    dummy.trend.fit <- list(trend = function(x) rep(NA_real_, length(x)),
                            std.dev = NA_real_)

    collected <- vector("list", ncol(x.means))
    for (i in seq_along(collected)) {
        fm <- fit.means[, i]
        # CV2 = variance / mean^2
        fcv2 <- fit.vars[, i] / fm^2
        if (ncells[i] >= 2L) {
            fit <- fitTrendCV2(fm, fcv2, ncells[i], ...)
        } else {
            fit <- dummy.trend.fit
        }

        xm <- unname(x.means[, i])
        # CV2 = variance / mean^2
        xcv2 <- unname(x.vars[, i]) / xm^2
        output <- DataFrame(mean = xm, total = xcv2, trend = fit$trend(xm))

        # Ratio of total CV2 to trend
        output$ratio <- output$total / output$trend
        # P-value: test if ratio > 1 (more variable than expected)
        output$p.value <- pnorm(output$ratio, mean = 1, sd = fit$std.dev,
                                lower.tail = FALSE)
        output$FDR <- p.adjust(output$p.value, method = "BH")

        rownames(output) <- rownames(x.means)
        metadata(output) <- c(list(mean = fm, cv2 = fcv2), fit)
        collected[[i]] <- output
    }
    names(collected) <- colnames(x.means)
    collected
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### scoreMarkers
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scran scoreMarkers computeMinRank
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom S4Vectors DataFrame SimpleList
setMethod("scoreMarkers", "DuckDBMatrix",
function(x, groups, block = NULL, pairings = NULL, lfc = 0, row.data = NULL,
         full.stats = FALSE, subset.row = NULL, BPPARAM = SerialParam(),
         true.auc = FALSE)
{
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Setup groups
    ncells <- ncol(x)
    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }
    groups <- as.factor(groups)
    group_levels <- levels(groups)
    ngroups <- nlevels(groups)

    # Subset rows if requested
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
        if (!is.null(row.data)) {
            row.data <- row.data[subset.row, , drop = FALSE]
        }
    }

    ngenes <- nrow(x)
    rnames <- rownames(x)

    # Compute per-group statistics via SQL (OPTIMIZED: combined query)
    # Previously this was 2 separate function calls (4 SQL queries)
    # Now it's 1 function call (2 SQL queries)
    stats <- .compute_group_stats_with_detection_DuckDBMatrix(x, groups,
                                                              block = NULL,
                                                              subset.row = NULL)

    # Determine pairings
    if (is.null(pairings)) {
        # All pairwise comparisons
        pairs <- expand.grid(left = seq_len(ngroups), right = seq_len(ngroups))
        pairs <- pairs[pairs$left != pairs$right, , drop = FALSE]
    } else {
        pairs <- .expand_pairings_DuckDBMatrix(pairings, group_levels)
    }

    # Compute effect sizes for each pair
    collapse.symmetric <- (lfc == 0)
    effect_sizes <- .compute_effect_sizes_DuckDBMatrix(
        stats$means, stats$vars, stats$ncells,
        stats$ndetected, stats$ncells,
        group_levels, pairs, lfc, collapse.symmetric
    )

    # Optionally compute true rank-based AUC via SQL window functions
    if (true.auc) {
        true_auc_matrix <- .compute_true_auc_sql_DuckDBMatrix(
            x, groups, group_levels, pairs, lfc
        )
        effect_sizes$auc <- true_auc_matrix
    }

    # Summarize effect sizes for each group
    result <- .summarize_effect_sizes_DuckDBMatrix(
        effect_sizes, group_levels, pairs, ngenes, rnames, full.stats
    )

    # Compute min-rank for each effect size
    for (g in names(result)) {
        rank_cols <- grep("^rank\\.", colnames(result[[g]]))
        for (j in rank_cols) {
            result[[g]][[j]] <- computeMinRank(result[[g]][[j]])
        }
    }

    # Add row data if provided
    if (!is.null(row.data)) {
        for (g in names(result)) {
            result[[g]] <- cbind(row.data, result[[g]])
        }
    }

    SimpleList(result)
})

# Compute true rank-based AUC via SQL window functions
# This uses the Wilcoxon rank-sum approach:
# AUC = (R1 - n1*(n1+1)/2) / (n1*n2)
# where R1 is the sum of ranks in group 1
#' @importFrom DBI dbGetQuery
#' @importFrom DuckDBDataFrame tblconn
.compute_true_auc_sql_DuckDBMatrix <-
function(x, groups, group_levels, pairs, lfc)
{
    # Get DuckDB connection and table info
    tbl <- tblconn(x)
    con <- dbplyr::remote_con(tbl)
    # Use sql_render to get the base SQL (dbplyr::remote_name returns NULL for complex queries)
    base_sql <- as.character(dbplyr::sql_render(tbl))

    # Extract table metadata using keynames for column name strings
    meta <- .extract_table_metadata(x)
    row_col <- meta$row_key
    col_col <- meta$col_key
    # The base_sql from dbplyr::sql_render() already handles data transformations
    # (e.g., normalizeCounts) and outputs the result as "value".
    # So we always use "value" as the data column name.
    data_col <- "value"

    ngenes <- nrow(x)
    npairs <- nrow(pairs)

    # Register groups as temporary table (zero-copy)
    group_df <- data.frame(
        cell_id = seq_along(groups),
        group_label = as.character(groups)
    )
    group_tbl_name <- paste0("auc_groups_", format(Sys.time(), "%H%M%S"),
                             "_", sample.int(10000, 1))
    duckdb::duckdb_register(con, group_tbl_name, group_df, overwrite = TRUE)
    on.exit(duckdb::duckdb_unregister(con, group_tbl_name), add = TRUE)

    # Initialize result matrix
    auc_matrix <- matrix(NA_real_, nrow = ngenes, ncol = npairs)

    # Compute AUC for each pairwise comparison
    for (p in seq_len(npairs)) {
        left_idx <- pairs$left[p]
        right_idx <- pairs$right[p]
        left_group <- group_levels[left_idx]
        right_group <- group_levels[right_idx]

        # Get cell counts for this pair
        n_left <- sum(groups == left_group)
        n_right <- sum(groups == right_group)

        # SQL for rank-based AUC computation
        # Uses window functions: O(n log n) per gene instead of O(n1*n2)
        lfc_expr <- if (lfc != 0) sprintf("- %f", lfc) else ""
        
        sql <- sprintf("
WITH data_tbl AS (%s),
all_cells_pair AS (
    SELECT cell_id, group_label as grp
    FROM %s
    WHERE group_label IN ('%s', '%s')
),
all_genes AS (
    SELECT DISTINCT %s as gene FROM data_tbl
),
-- Full outer join to get all (gene, cell) combinations with zeros for sparse data
full_data AS (
    SELECT
        g.gene,
        c.cell_id as cell,
        COALESCE(d.%s, 0) %s as val,
        c.grp
    FROM all_genes g
    CROSS JOIN all_cells_pair c
    LEFT JOIN data_tbl d ON g.gene = d.%s AND c.cell_id = d.%s
),
-- Assign row numbers for ranking (ordered by value, then cell for tie-breaking)
with_rownum AS (
    SELECT
        gene, cell, val, grp,
        ROW_NUMBER() OVER (PARTITION BY gene ORDER BY val, cell) as row_num
    FROM full_data
),
-- Compute average rank for tied values within each gene
with_avg_rank AS (
    SELECT
        gene, cell, val, grp,
        AVG(row_num) OVER (PARTITION BY gene, val) as avg_rank
    FROM with_rownum
),
-- Sum ranks for the left group
rank_sums AS (
    SELECT
        gene,
        SUM(CASE WHEN grp = '%s' THEN avg_rank ELSE 0 END) as R1
    FROM with_avg_rank
    GROUP BY gene
)
SELECT
    gene,
    (R1 - %d * (%d + 1.0) / 2.0) / (%d * %d) as auc
FROM rank_sums
ORDER BY gene
",
            base_sql,
            group_tbl_name,
            left_group, right_group,
            row_col,
            data_col, lfc_expr,
            row_col, col_col,
            left_group,
            n_left, n_left, n_left, n_right
        )

        result <- DBI::dbGetQuery(con, sql)

        # Map results back to matrix (handle potentially missing genes)
        gene_idx <- match(result$gene, seq_len(ngenes))
        if (any(!is.na(gene_idx))) {
            auc_matrix[gene_idx[!is.na(gene_idx)], p] <- result$auc[!is.na(gene_idx)]
        }
    }

    auc_matrix
}

.expand_pairings_DuckDBMatrix <- function(pairings, group_levels)
{
    if (is.vector(pairings) || is.factor(pairings)) {
        # Subset of groups for all pairwise comparisons
        idx <- match(as.character(pairings), group_levels)
        pairs <- expand.grid(left = idx, right = idx)
        pairs <- pairs[pairs$left != pairs$right, , drop = FALSE]
    } else if (is.list(pairings) && length(pairings) == 2L) {
        # Two vectors defining left and right groups
        left_idx <- match(as.character(pairings[[1L]]), group_levels)
        right_idx <- match(as.character(pairings[[2L]]), group_levels)
        pairs <- expand.grid(left = left_idx, right = right_idx)
        pairs <- pairs[pairs$left != pairs$right, , drop = FALSE]
    } else if (is.matrix(pairings)) {
        # Explicit pairs
        pairs <- data.frame(
            left = match(pairings[, 1L], group_levels),
            right = match(pairings[, 2L], group_levels)
        )
    } else {
        stop("invalid 'pairings' specification")
    }
    pairs
}

.compute_effect_sizes_DuckDBMatrix <-
function(means, vars, ncells, ndetected, det_ncells, group_levels, pairs,
         lfc, collapse.symmetric)
{
    # VECTORIZED implementation - no for-loops
    # Extract left and right group names for all pairs at once
    left_names <- group_levels[pairs$left]
    right_names <- group_levels[pairs$right]

    # Extract matrices of left/right statistics using column indexing
    # Each column corresponds to a pair, each row to a gene
    left_means <- means[, left_names, drop = FALSE]
    right_means <- means[, right_names, drop = FALSE]
    left_vars <- vars[, left_names, drop = FALSE]
    right_vars <- vars[, right_names, drop = FALSE]

    # Cell counts per group (vector, one per pair)
    left_n <- ncells[left_names]
    right_n <- ncells[right_names]

    # Detection counts (matrix: genes x pairs)
    left_det <- ndetected[, left_names, drop = FALSE]
    right_det <- ndetected[, right_names, drop = FALSE]
    left_det_n <- det_ncells[left_names]
    right_det_n <- det_ncells[right_names]

    # Log-fold change (delta mean) - vectorized
    delta_mean <- left_means - right_means - lfc

    # Cohen's d (standardized mean difference) - vectorized
    pooled_var <- (left_vars + right_vars) / 2
    pooled_var <- pmax(pooled_var, 1e-8)
    cohen_d <- (left_means - right_means - lfc) / sqrt(pooled_var)

    # AUC (probability that left > right) - vectorized
    # NOTE: This uses a normal approximation, not true rank-based AUC.
    # For true rank-based statistics, use findMarkers(test.type="wilcox")
    # Broadcast cell counts across rows (genes)
    # left_n and right_n are vectors of length npairs
    combined_var <- sweep(left_vars, 2, left_n, `/`) +
                    sweep(right_vars, 2, right_n, `/`)
    combined_var <- pmax(combined_var, 1e-8)
    z <- (left_means - right_means - lfc) / sqrt(combined_var)
    auc <- pnorm(z)

    # Log-fold change in detection proportion - vectorized
    # Broadcast detection cell counts across rows (genes)
    left_prop <- sweep(left_det + 0.5, 2, left_det_n + 1, `/`)
    right_prop <- sweep(right_det + 0.5, 2, right_det_n + 1, `/`)
    lfc_detected <- log2(left_prop) - log2(right_prop)

    # Convert to plain matrices (drop dimnames for consistency)
    list(
        cohen_d = unname(as.matrix(cohen_d)),
        auc = unname(as.matrix(auc)),
        delta_mean = unname(as.matrix(delta_mean)),
        lfc_detected = unname(as.matrix(lfc_detected)),
        pairs = pairs
    )
}

#' @importFrom S4Vectors DataFrame
#' @importFrom stats median
.summarize_effect_sizes_DuckDBMatrix <-
function(effect_sizes, group_levels, pairs, ngenes, rnames, full.stats)
{
    ngroups <- length(group_levels)
    result <- list()

    for (g in seq_len(ngroups)) {
        group_name <- group_levels[g]

        # Find pairs involving this group on the left
        left_idx <- which(pairs$left == g)

        if (length(left_idx) == 0L) {
            # No comparisons for this group
            df <- DataFrame(
                self.average = rep(NA_real_, ngenes),
                other.average = rep(NA_real_, ngenes),
                self.detected = rep(NA_real_, ngenes),
                other.detected = rep(NA_real_, ngenes)
            )
            rownames(df) <- rnames
            result[[group_name]] <- df
            next
        }

        # Extract effect sizes for this group's comparisons
        cohen_subset <- effect_sizes$cohen_d[, left_idx, drop = FALSE]
        auc_subset <- effect_sizes$auc[, left_idx, drop = FALSE]
        delta_subset <- effect_sizes$delta_mean[, left_idx, drop = FALSE]
        det_subset <- effect_sizes$lfc_detected[, left_idx, drop = FALSE]

        # Compute summary statistics using matrixStats for performance
        # (rowMedians, rowMins, rowMaxs are much faster than apply())
        df <- DataFrame(
            # Cohen's d summaries
            mean.logFC.cohen = rowMeans(cohen_subset, na.rm = TRUE),
            median.logFC.cohen = matrixStats::rowMedians(cohen_subset, na.rm = TRUE),
            min.logFC.cohen = matrixStats::rowMins(cohen_subset, na.rm = TRUE),
            max.logFC.cohen = matrixStats::rowMaxs(cohen_subset, na.rm = TRUE),
            rank.logFC.cohen = cohen_subset,  # For computeMinRank later

            # AUC summaries
            mean.AUC = rowMeans(auc_subset, na.rm = TRUE),
            median.AUC = matrixStats::rowMedians(auc_subset, na.rm = TRUE),
            min.AUC = matrixStats::rowMins(auc_subset, na.rm = TRUE),
            max.AUC = matrixStats::rowMaxs(auc_subset, na.rm = TRUE),
            rank.AUC = auc_subset,

            # Delta mean summaries
            mean.logFC.detected = rowMeans(det_subset, na.rm = TRUE),
            median.logFC.detected = matrixStats::rowMedians(det_subset, na.rm = TRUE),
            min.logFC.detected = matrixStats::rowMins(det_subset, na.rm = TRUE),
            max.logFC.detected = matrixStats::rowMaxs(det_subset, na.rm = TRUE),
            rank.logFC.detected = det_subset
        )

        if (full.stats) {
            # Add full pairwise statistics
            other_groups <- group_levels[pairs$right[left_idx]]
            colnames(cohen_subset) <- paste0("full.logFC.cohen.", other_groups)
            colnames(auc_subset) <- paste0("full.AUC.", other_groups)
            colnames(det_subset) <- paste0("full.logFC.detected.", other_groups)

            df <- cbind(df,
                        as(cohen_subset, "DataFrame"),
                        as(auc_subset, "DataFrame"),
                        as(det_subset, "DataFrame"))
        }

        rownames(df) <- rnames
        result[[group_name]] <- df
    }

    result
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summaryMarkerStats
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scran summaryMarkerStats
#' @importFrom scuttle .bpNotSharedOrUp
#' @importFrom S4Vectors DataFrame SimpleList
setMethod("summaryMarkerStats", "DuckDBMatrix",
function(x, groups, row.data = NULL, average = "mean", BPPARAM = SerialParam())
{
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    # Setup groups
    ncells <- ncol(x)
    if (length(groups) != ncells) {
        stop("length of 'groups' does not equal 'ncol(x)'")
    }
    groups <- as.factor(groups)
    group_levels <- levels(groups)
    ngroups <- nlevels(groups)

    ngenes <- nrow(x)
    rnames <- rownames(x)

    # Compute per-group means AND detection proportions via SQL (OPTIMIZED)
    # Previously this was 2 separate function calls (4 SQL queries)
    # Now it's 1 function call (2 SQL queries)
    stats <- .compute_group_stats_with_detection_DuckDBMatrix(x, groups,
                                                              block = NULL,
                                                              subset.row = NULL)

    # Build per-group result DataFrames
    result <- list()
    for (g in seq_len(ngroups)) {
        group_name <- group_levels[g]

        # Self statistics
        self_avg <- stats$means[, group_name]
        self_det <- stats$ndetected[, group_name] / stats$ncells[group_name]

        # Other statistics (mean across other groups)
        other_groups <- group_levels[-g]
        if (length(other_groups) > 0L) {
            other_avg <- rowMeans(stats$means[, other_groups, drop = FALSE])
            other_det <- rowMeans(
                sweep(stats$ndetected[, other_groups, drop = FALSE],
                      2, stats$ncells[other_groups], "/")
            )
        } else {
            other_avg <- rep(NA_real_, ngenes)
            other_det <- rep(NA_real_, ngenes)
        }

        df <- DataFrame(
            self.average = self_avg,
            other.average = other_avg,
            self.detected = self_det,
            other.detected = other_det
        )
        rownames(df) <- rnames

        # Add row data if provided
        if (!is.null(row.data)) {
            df <- cbind(row.data, df)
        }

        result[[group_name]] <- df
    }

    SimpleList(result)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### correlatePairs
###

#' @export
#' @importClassesFrom DuckDBArray DuckDBMatrix
#' @importFrom BiocParallel SerialParam
#' @importFrom scran correlatePairs
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr collect group_by summarize
#' @importMethodsFrom scran correlatePairs
setMethod("correlatePairs", "DuckDBMatrix",
function(x, subset.row = NULL, pairings = NULL, use.names = TRUE,
         BPPARAM = SerialParam(), ...)
{
    # Only support fill = 0 (standard sparse representation)
    fill <- x@seed@fill
    if (fill != 0) {
        stop("correlatePairs for DuckDBMatrix only supports fill = 0")
    }

    # Require subset.row for large matrices
    ngenes <- nrow(x)
    if (is.null(subset.row) && ngenes > 1000L) {
        stop("subset.row is required when nrow(x) > 1000 to avoid ",
             "computing too many pairs. Use subset.row to select genes ",
             "of interest (e.g., highly variable genes).")
    }

    # Subset if requested
    if (!is.null(subset.row)) {
        x <- x[subset.row, , drop = FALSE]
    }

    ngenes <- nrow(x)
    ncells <- ncol(x)
    rnames <- rownames(x)

    if (ngenes < 2L) {
        stop("need at least 2 genes to compute pairwise correlations")
    }

    # Extract table metadata (REFACTORED)
    meta <- .extract_table_metadata(x)
    table <- meta$table
    row_key <- meta$row_key
    col_key <- meta$col_key
    row_keycol <- meta$row_keycol
    col_keycol <- meta$col_keycol
    datacol <- meta$datacol

    # Compute row statistics: sum, sum_sq for each gene
    row_stats <- .compute_row_stats_for_cor(table, row_key, datacol)

    # Compute pairwise sum of products via self-join
    pair_stats <- .compute_pair_products(table, row_key, col_key, datacol,
                                         pairings, row_keycol)

    # Compute correlations using sparse-aware formula
    result <- .compute_sparse_pearson(row_stats, pair_stats, row_key,
                                       row_keycol, ncells, pairings)

    # Format output
    if (use.names && !is.null(rnames)) {
        result$gene1 <- rnames[match(result$gene1, row_keycol)]
        result$gene2 <- rnames[match(result$gene2, row_keycol)]
    }

    DataFrame(result)
})

#' @importFrom BiocGenerics as.data.frame
#' @importFrom dplyr collect group_by summarize
#' @importFrom DuckDBDataFrame tblconn
.compute_row_stats_for_cor <- function(table, row_key, datacol)
{
    # Use select = TRUE to get the transformed datacol values
    conn <- tblconn(table)
    row_sym <- as.name(row_key)
    datacol_name <- as.name(names(table@datacols)[1L])

    # Compute sum and sum of squares per gene
    aggr <- list(
        sum_val = call("sum", datacol_name, na.rm = TRUE),
        sum_sq = call("sum", call("*", datacol_name, datacol_name), na.rm = TRUE)
    )

    result <- conn
    result <- group_by(result, !!row_sym)
    result <- summarize(result, !!!aggr, .groups = "drop")
    result <- collect(result)

    as.data.frame(result)
}

#' @importFrom BiocGenerics as.data.frame
#' @importFrom dplyr collect group_by summarize sql
#' @importFrom DuckDBDataFrame tblconn
.compute_pair_products <- function(table, row_key, col_key, datacol,
                                    pairings, row_keycol)
{
    # Use select = TRUE to get the transformed datacol values
    conn <- tblconn(table)

    # Build self-join query for pairwise products
    # SELECT a.gene as gene1, b.gene as gene2, SUM(a.value * b.value) as sum_prod
    # FROM data a JOIN data b ON a.cell = b.cell AND a.gene < b.gene
    # GROUP BY a.gene, b.gene

    # Get underlying connection
    db_conn <- conn$src$con

    # Build SQL directly for self-join (dbplyr doesn't handle self-joins well)
    # Get the SQL for the base table
    base_sql <- dbplyr::sql_render(conn)

    if (is.null(pairings)) {
        # All pairs
        query <- sprintf("
            SELECT a.\"%s\" as gene1, b.\"%s\" as gene2,
                   SUM(a.\"%s\" * b.\"%s\") as sum_prod
            FROM (%s) a
            JOIN (%s) b ON a.\"%s\" = b.\"%s\" AND a.\"%s\" < b.\"%s\"
            GROUP BY a.\"%s\", b.\"%s\"
        ", row_key, row_key,
           names(table@datacols)[1L], names(table@datacols)[1L],
           base_sql, base_sql,
           col_key, col_key, row_key, row_key,
           row_key, row_key)
    } else {
        # Specific pairs - filter to requested pairs
        # Convert pairings to gene keys if needed
        if (is.character(pairings)) {
            pairings <- matrix(match(pairings, names(row_keycol)), ncol = 2)
        }
        pairings <- matrix(row_keycol[pairings], ncol = 2)

        # Create pairs table
        pairs_values <- paste(sprintf("(%s, %s)", pairings[, 1], pairings[, 2]),
                              collapse = ", ")

        query <- sprintf("
            WITH pairs AS (SELECT * FROM (VALUES %s) AS t(g1, g2))
            SELECT a.\"%s\" as gene1, b.\"%s\" as gene2,
                   SUM(a.\"%s\" * b.\"%s\") as sum_prod
            FROM (%s) a
            JOIN (%s) b ON a.\"%s\" = b.\"%s\"
            JOIN pairs p ON (a.\"%s\" = p.g1 AND b.\"%s\" = p.g2)
                         OR (a.\"%s\" = p.g2 AND b.\"%s\" = p.g1)
            WHERE a.\"%s\" < b.\"%s\"
            GROUP BY a.\"%s\", b.\"%s\"
        ", pairs_values,
           row_key, row_key,
           names(table@datacols)[1L], names(table@datacols)[1L],
           base_sql, base_sql,
           col_key, col_key,
           row_key, row_key, row_key, row_key,
           row_key, row_key,
           row_key, row_key)
    }

    result <- DBI::dbGetQuery(db_conn, query)
    as.data.frame(result)
}

.compute_sparse_pearson <- function(row_stats, pair_stats, row_key,
                                     row_keycol, ncells, pairings)
{
    k <- as.numeric(ncells)

    # Create lookup for row stats
    row_stats_lookup <- row_stats
    names(row_stats_lookup)[names(row_stats_lookup) == row_key] <- "gene"

    # If no pairs found (all genes have disjoint non-zeros), create all pairs
    if (nrow(pair_stats) == 0L) {
        if (is.null(pairings)) {
            # Generate all pairs
            genes <- row_keycol
            n <- length(genes)
            pairs <- expand.grid(gene1 = genes, gene2 = genes)
            pairs <- pairs[pairs$gene1 < pairs$gene2, , drop = FALSE]
            pair_stats <- data.frame(gene1 = pairs$gene1, gene2 = pairs$gene2,
                                     sum_prod = 0)
        } else {
            # Use specified pairings with sum_prod = 0
            pair_stats <- data.frame(gene1 = pairings[, 1], gene2 = pairings[, 2],
                                     sum_prod = 0)
        }
    }

    # Handle missing pairs (genes with no overlapping non-zeros)
    if (is.null(pairings)) {
        genes <- unique(c(pair_stats$gene1, pair_stats$gene2))
        all_genes <- row_keycol
        genes <- all_genes[all_genes %in% row_stats_lookup$gene]
        n <- length(genes)

        # Check if all pairs are present
        expected_pairs <- n * (n - 1L) / 2L
        if (nrow(pair_stats) < expected_pairs) {
            # Add missing pairs with sum_prod = 0
            existing <- paste(pair_stats$gene1, pair_stats$gene2, sep = "_")
            all_pairs <- expand.grid(gene1 = genes, gene2 = genes)
            all_pairs <- all_pairs[all_pairs$gene1 < all_pairs$gene2, , drop = FALSE]
            all_pairs$key <- paste(all_pairs$gene1, all_pairs$gene2, sep = "_")
            missing <- all_pairs[!all_pairs$key %in% existing, c("gene1", "gene2")]
            if (nrow(missing) > 0L) {
                missing$sum_prod <- 0
                pair_stats <- rbind(pair_stats, missing)
            }
        }
    }

    # Merge row stats for gene1 and gene2
    pair_stats <- merge(pair_stats,
                        setNames(row_stats_lookup, c("gene1", "sum1", "sum_sq1")),
                        by = "gene1", all.x = TRUE)
    pair_stats <- merge(pair_stats,
                        setNames(row_stats_lookup, c("gene2", "sum2", "sum_sq2")),
                        by = "gene2", all.x = TRUE)

    # Compute Pearson correlation using sparse formula:
    # r = (sum_prod - sum1 * sum2 / k) / sqrt(SS1 * SS2)
    # where SS = sum_sq - sum^2 / k
    pair_stats$SS1 <- pair_stats$sum_sq1 - pair_stats$sum1^2 / k
    pair_stats$SS2 <- pair_stats$sum_sq2 - pair_stats$sum2^2 / k

    numerator <- pair_stats$sum_prod - pair_stats$sum1 * pair_stats$sum2 / k
    denominator <- sqrt(pair_stats$SS1 * pair_stats$SS2)

    # Handle edge cases (zero variance)
    pair_stats$rho <- ifelse(denominator > 0, numerator / denominator, NA_real_)

    # Return ordered by gene1, gene2
    pair_stats <- pair_stats[order(pair_stats$gene1, pair_stats$gene2),
                             c("gene1", "gene2", "rho"), drop = FALSE]
    rownames(pair_stats) <- NULL

    pair_stats
}
