#' DuckDBMatrix scran methods
#'
#' @description
#' scran methods for \link[DuckDBArray:DuckDBMatrix-class]{DuckDBMatrix} objects.
#'
#' @section Pairwise Correlation Methods:
#' The following pairwise correlation methods have optimized DuckDBMatrix implementations:
#' \describe{
#'   \item{\code{correlatePairs(x, subset.row = NULL, pairings = NULL,
#'     use.names = TRUE, BPPARAM = SerialParam())}:}{
#'     Compute pairwise Pearson correlations between genes using sparse-aware
#'     SQL aggregation, with a standard t-test significance for each
#'     correlation (\code{p.value}, BH-adjusted as \code{FDR}). Only
#'     \code{fill = 0} is supported. This is Pearson's r on \code{x} as given
#'     (typically log-normalized values), not \pkg{scran}'s Spearman's rho on
#'     ranked expression -- results are not expected to numerically match
#'     \code{scran::correlatePairs()}.
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
#'   \item \code{\link[DuckDBArray:DuckDBMatrix-class]{DuckDBMatrix-class}} for the main class
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
#' @return
#' These methods mirror their \pkg{scran} generics. \code{correlatePairs()}
#' returns a \code{DataFrame} of gene pairs with columns \code{gene1}, \code{gene2},
#' \code{rho} (Pearson's r), \code{p.value}, and \code{FDR}.
#' \code{modelGeneVar()}, \code{modelGeneVarByPoisson()}, and
#' \code{modelGeneCV2()} return a \code{DataFrame} of per-gene variance
#' statistics, with the fitted mean-variance trend in \code{metadata()}.
#' \code{pairwiseTTests()} and \code{pairwiseBinom()} return a list with
#' \code{statistics} and \code{pairs} elements; \code{findMarkers()} and
#' \code{scoreMarkers()} return a \link[S4Vectors]{List} of per-group
#' \code{DataFrame}s of marker statistics; and \code{summaryMarkerStats()}
#' returns a \code{DataFrame} summarizing marker statistics per group.
#'
#' @aliases correlatePairs,DuckDBMatrix-method
#' @aliases modelGeneVar,DuckDBMatrix-method
#' @aliases modelGeneVarByPoisson,DuckDBMatrix-method
#' @aliases modelGeneCV2,DuckDBMatrix-method
#' @aliases pairwiseTTests,DuckDBMatrix-method
#' @aliases pairwiseBinom,DuckDBMatrix-method
#' @aliases findMarkers,DuckDBMatrix-method
#' @aliases scoreMarkers,DuckDBMatrix-method
#' @aliases summaryMarkerStats,DuckDBMatrix-method
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
#' head(scran::modelGeneVarByPoisson(pqmat))
#' unlink(tmp, recursive = TRUE)
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
.temp_table_counter_env <- new.env(parent = emptyenv())
.temp_table_counter_env$n <- 0L

.next_temp_table_suffix <- function() {
    .temp_table_counter_env$n <- .temp_table_counter_env$n + 1L
    paste0(Sys.getpid(), "_", .temp_table_counter_env$n)
}

.register_temp_table <- function(db_conn, df, prefix = "temp") {
    table_name <- paste0(prefix, "_", .next_temp_table_suffix())
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

# Shared primitive: fused single-pass computation of per-gene means and
# variances within named groups (blocks, or group x block combinations),
# including the implicit "fill" population (cells absent from the sparse
# table).
#
# This used to be a numerically-stable but wasteful two-pass algorithm
# (Pass 1: SQL to get per-group means; Pass 2: re-join those means and sum
# squared deviations in SQL). It is now fused into ONE GROUP BY query that
# computes, per (gene, group), only quantities over the STORED rows:
#   n = COUNT(*), sum_y = SUM(y), S2 = COALESCE(VAR_SAMP(y), 0)
# using DuckDB's own numerically stable VAR_SAMP -- never a hand-rolled
# sum-of-squares. The implicit "fill" population (cells absent from the
# sparse table, k - n of them per group, where k is the group size known
# from the assignment rather than the database) is then folded in
# afterwards, in R, via the Chan et al. parallel-variance-combination
# identity:
#
#   Sigma(y - mu)^2 over k  =  Sigma(y - ybar)^2 over n
#                            + n * (ybar - mu)^2
#                            + (k - n) * (fill - mu)^2
#
# where mu is the TRUE overall mean including the implicit fill cells:
#   mu = (sum_y + fill * (k - n)) / k
# The stable (n - 1) * S2 term comes straight from DuckDB's VAR_SAMP; the two
# correction terms are each a small count times a squared difference of
# already-stable means, so no raw sum-of-squares of individual values is
# ever computed. Setting fill = 0 recovers exactly the single-pass fast path
# already used by DuckDBArray's `.marginVars()` fill == 0 branch (see
# DuckDBTable-matrixStats.R) -- this is the same algebra, generalized to a
# nonzero fill.
#
# `skip_empty_groups` controls behaviour for a group with k == 0 cells: when
# FALSE (blocks, which by construction cover every cell), the division by k
# is allowed to run and naturally yields NaN, matching the historical
# blocked-stats behaviour for an empty block level; when TRUE (groups, which
# can be empty after NA exclusion), the group is skipped and its mean/
# variance stay NA, matching the historical group-stats behaviour.
#' @importFrom dplyr collect group_by inner_join left_join mutate summarize tbl
.fused_group_variance_stats <-
function(conn, db_conn, row_key, col_key, value_expr, assign_df, id_col,
         col_ids, group_sizes, col_names, ngenes, rnames, row_keyvals,
         fill, inner = FALSE, mutate_expr = NULL, skip_empty_groups = FALSE)
{
    assign_tbl_name <- .register_temp_table(db_conn, assign_df, "assign")
    on.exit(.unregister_temp_tables(db_conn, assign_tbl_name), add = TRUE)

    row_sym <- as.name(row_key)
    id_sym <- as.name(id_col)

    joined <- if (inner) {
        inner_join(conn, tbl(db_conn, assign_tbl_name), by = col_key)
    } else {
        left_join(conn, tbl(db_conn, assign_tbl_name), by = col_key)
    }
    if (!is.null(mutate_expr)) {
        joined <- mutate(joined, !!!mutate_expr)
    }

    # Single fused GROUP BY: n, sum, and DuckDB's native VAR_SAMP together.
    aggr <- list(
        stored_n = call("n"),
        sum_y = call("sum", value_expr, na.rm = TRUE),
        S2 = call("coalesce", call("var_samp", value_expr), 0)
    )
    result <- group_by(joined, !!row_sym, !!id_sym)
    result <- summarize(result, !!!aggr, .groups = "drop")
    result <- collect(result)

    # GOTCHA: DuckDB returns COUNT(*)/n() as a bit64::integer64 via collect().
    # Arithmetic on the uncast integer64 silently yields WRONG variances (not
    # float noise -- ~1.0-1.4 absolute error was measured), so cast to
    # double immediately, before any arithmetic below.
    result$stored_n <- as.double(result$stored_n)

    means <- matrix(NA_real_, nrow = ngenes, ncol = length(col_ids),
                    dimnames = list(rnames, col_names))
    vars <- matrix(NA_real_, nrow = ngenes, ncol = length(col_ids),
                   dimnames = list(rnames, col_names))

    for (pos in seq_along(col_ids)) {
        g <- col_ids[pos]
        k <- group_sizes[[pos]]
        if (skip_empty_groups && k == 0L) next

        sub <- result[result[[id_col]] == g, , drop = FALSE]
        n <- numeric(ngenes)
        sum_y <- numeric(ngenes)
        S2 <- numeric(ngenes)
        if (nrow(sub) > 0L) {
            row_idx <- match(sub[[row_key]], row_keyvals)
            n[row_idx] <- sub$stored_n
            sum_y[row_idx] <- sub$sum_y
            S2[row_idx] <- sub$S2
        }

        mu <- (sum_y + fill * (k - n)) / k
        means[, pos] <- mu

        if (k <= 1L) {
            vars[, pos] <- NA_real_
            next
        }

        ybar <- ifelse(n > 0, sum_y / n, fill)
        sum_sq_dev_full <- (n - 1) * S2 + n * (ybar - mu)^2 + (k - n) * (fill - mu)^2
        vars[, pos] <- sum_sq_dev_full / (k - 1)
    }

    list(means = means, vars = vars)
}

# Compute per-gene means and variances, optionally by block.
# Backs modelGeneVar()/modelGeneVarByPoisson(). Delegates to the shared
# fused single-pass primitive, `.fused_group_variance_stats()` (see above).
#' @importFrom DuckDBDataFrame tblconn
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
    row_keyvals <- unname(row_keycol)

    stats <- .fused_group_variance_stats(
        conn, db_conn, row_key, col_key, datacol, block_df, "block_id",
        col_ids = levels(block), group_sizes = block_ncells,
        col_names = bnames, ngenes = ngenes, rnames = rnames,
        row_keyvals = row_keyvals, fill = fill,
        inner = FALSE, skip_empty_groups = FALSE
    )

    list(means = stats$means, vars = stats$vars, ncells = block_ncells)
}

# Compute per-gene means and variances by group (for t-tests).
# Backs pairwiseTTests()/findMarkers(). Delegates to the shared fused
# single-pass primitive, `.fused_group_variance_stats()` (see above
# `.compute_blocked_stats_DuckDBMatrix` for the algorithm derivation).
#' @importFrom DuckDBDataFrame tblconn
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
    row_keyvals <- unname(row_keycol)

    # group_df already contains exactly the cells with a non-NA group
    # assignment, so an inner join on it is equivalent to the historical
    # filter(col_key %in% group_df[[col_key]]) |> left_join(...): no row can
    # pass that filter without also matching group_df, so the inner join
    # subsumes it (and avoids ever materializing a large literal IN (...)
    # list in the generated SQL).
    stats <- .fused_group_variance_stats(
        conn, db_conn, row_key, col_key, datacol, group_df, "group_id",
        col_ids = combined_levels, group_sizes = group_sizes,
        col_names = combined_levels, ngenes = ngenes, rnames = rnames,
        row_keyvals = row_keyvals, fill = fill,
        inner = TRUE, skip_empty_groups = TRUE
    )

    if (is.null(block)) {
        list(means = stats$means, vars = stats$vars, ncells = group_sizes)
    } else {
        # Reshape into block structure
        list(means = stats$means, vars = stats$vars, ncells = group_sizes,
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
                warning("no within-block comparison between ", host,
                        " and ", target)
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

    # Simulate Poisson counts and compute log-normalized variance. All npts x
    # ncells draws happen in one vectorized call rather than one rpois()/
    # rnbinom() call per simulation point.
    ncells <- length(size.factors)
    lambda_mat <- outer(sim_means, size.factors)
    if (dispersion > 0) {
        # Negative binomial
        counts_mat <- matrix(rnbinom(length(lambda_mat), mu = lambda_mat,
                                     size = 1 / dispersion),
                             nrow = npts, ncol = ncells)
    } else {
        # Poisson
        counts_mat <- matrix(rpois(length(lambda_mat), lambda = lambda_mat),
                             nrow = npts, ncol = ncells)
    }
    # Log-normalize (each column divided by its own cell's size factor)
    log_vals_mat <- log2(sweep(counts_mat, 2L, size.factors, "/") + pseudo.count)
    sim_vars <- matrixStats::rowVars(log_vals_mat)

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
# (operates on the derived column value/sf). Backs modelGeneCV2(). Delegates
# to the shared fused single-pass primitive, `.fused_group_variance_stats()`
# (see above `.compute_blocked_stats_DuckDBMatrix` for the algorithm
# derivation); the normalization division is folded in as a `mutate_expr`
# applied inside the same single GROUP BY query, so there is still only one
# physical table scan.
#' @importFrom DuckDBDataFrame tblconn
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
    row_keyvals <- unname(row_keycol)

    # Normalized fill value for the implicit "fill" population: since each
    # missing cell's true size factor is unknown, approximate with the mean
    # size factor across all cells (matches the historical implementation).
    norm_fill <- fill / mean(size.factors)

    # Compute normalized value (value / sf) as part of the single fused query.
    mutate_expr <- list(norm_val = call("/", datacol, as.name("sf")))

    stats <- .fused_group_variance_stats(
        conn, db_conn, row_key, col_key, as.name("norm_val"), sf_df, "block_id",
        col_ids = levels(block), group_sizes = block_ncells,
        col_names = bnames, ngenes = ngenes, rnames = rnames,
        row_keyvals = row_keyvals, fill = norm_fill,
        inner = FALSE, mutate_expr = mutate_expr, skip_empty_groups = FALSE
    )

    list(means = stats$means, vars = stats$vars, ncells = block_ncells)
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
    row_keyvals <- unname(meta$row_keycol)

    # Register groups as temporary table (zero-copy)
    group_df <- data.frame(
        cell_id = seq_along(groups),
        group_label = as.character(groups)
    )
    group_tbl_name <- .register_temp_table(con, group_df, "auc_groups")
    temp_tables <- group_tbl_name

    # Register the full universe of raw gene keys as a temporary table
    genes_df <- data.frame(gene = row_keyvals)
    genes_tbl_name <- .register_temp_table(con, genes_df, "auc_genes")
    temp_tables <- c(temp_tables, genes_tbl_name)
    on.exit(.unregister_temp_tables(con, temp_tables), add = TRUE)

    # Initialize result matrix
    auc_matrix <- matrix(NA_real_, nrow = ngenes, ncol = npairs)

    # Compute AUC for each pairwise comparison
    for (p in seq_len(npairs)) {
        left_idx <- pairs$left[p]
        right_idx <- pairs$right[p]
        left_group <- group_levels[left_idx]
        right_group <- group_levels[right_idx]
        # group labels are arbitrary user-supplied text (factor levels), so they
        # must be quoted/escaped rather than interpolated raw into the SQL below
        left_group_sql <- DBI::dbQuoteString(con, left_group)
        right_group_sql <- DBI::dbQuoteString(con, right_group)

        # Get cell counts for this pair
        n_left <- sum(groups == left_group)
        n_right <- sum(groups == right_group)

        # SQL for rank-based AUC computation
        # Uses window functions: O(n log n) per gene instead of O(n1*n2)
        # A threshold shift must move the two groups relative to each other
        # (P(left > right + lfc)), not shift every value by the same amount
        val_expr <- if (lfc != 0) {
            sprintf("CASE WHEN c.grp = %s THEN COALESCE(d.%s, 0) - %.10f ELSE COALESCE(d.%s, 0) END",
                    left_group_sql, data_col, lfc, data_col)
        } else {
            sprintf("COALESCE(d.%s, 0)", data_col)
        }

        sql <- sprintf("
WITH data_tbl AS (%s),
all_cells_pair AS (
    SELECT cell_id, group_label as grp
    FROM %s
    WHERE group_label IN (%s, %s)
),
all_genes AS (
    SELECT gene FROM %s
),
-- Full outer join to get all (gene, cell) combinations with zeros for sparse data
full_data AS (
    SELECT
        g.gene,
        c.cell_id as cell,
        %s as val,
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
        SUM(CASE WHEN grp = %s THEN avg_rank ELSE 0 END) as R1
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
            left_group_sql, right_group_sql,
            genes_tbl_name,
            val_expr,
            row_col, col_col,
            left_group_sql,
            n_left, n_left, n_left, n_right
        )

        result <- DBI::dbGetQuery(con, sql)

        # Map results back to matrix (handle potentially missing genes)
        gene_idx <- match(result$gene, row_keyvals)
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

    # Compute row statistics: n, sum, and DuckDB's native VAR_SAMP for each
    # gene, over the explicitly-stored cells only.
    row_stats <- .compute_row_stats_for_cor(table, row_key, datacol, row_keycol)

    # Compute pairwise "both explicitly stored" statistics via self-join:
    # n11, the two sub-sums, and DuckDB's native COVAR_SAMP.
    pair_stats <- .compute_pair_products(table, row_key, col_key, datacol,
                                         pairings, row_keycol)

    # Compute correlations using a numerically-stable sparse-aware formula
    result <- .compute_sparse_pearson(row_stats, pair_stats,
                                       row_keycol, ncells, pairings, fill)

    # Format output
    if (use.names && !is.null(rnames)) {
        result$gene1 <- rnames[match(result$gene1, row_keycol)]
        result$gene2 <- rnames[match(result$gene2, row_keycol)]
    }

    DataFrame(result)
})

#' @importFrom dplyr collect group_by summarize
#' @importFrom DuckDBDataFrame tblconn
.compute_row_stats_for_cor <- function(table, row_key, datacol, row_keycol)
{
    # Use select = TRUE to get the transformed datacol values
    conn <- tblconn(table)
    row_sym <- as.name(row_key)
    datacol_name <- as.name(names(table@datacols)[1L])

    # Per-gene statistics over the EXPLICITLY STORED cells only: n (count),
    # sum_val (SUM), and S2 (DuckDB's native, numerically stable VAR_SAMP).
    # Deliberately NOT SUM(value * value): see the long comment above
    # .compute_sparse_pearson() for why that raw-moment formula is unsafe.
    aggr <- list(
        n = call("n"),
        sum_val = call("sum", datacol_name, na.rm = TRUE),
        S2 = call("coalesce", call("var_samp", datacol_name), 0)
    )

    result <- conn
    result <- group_by(result, !!row_sym)
    result <- summarize(result, !!!aggr, .groups = "drop")
    result <- collect(result)
    result <- as.data.frame(result)

    # GOTCHA (same as .fused_group_variance_stats()): DuckDB returns
    # COUNT(*)/n() as bit64::integer64 via collect(). Cast to double before
    # any arithmetic, or downstream arithmetic silently goes wrong.
    result$n <- as.double(result$n)
    names(result)[names(result) == row_key] <- "gene"

    # A gene with ZERO explicitly-stored cells (an all-fill/all-zero row)
    # never produces a group in a GROUP BY result. Reindex over the full
    # requested gene set so every gene has a row (n = 0, sum_val = 0,
    # S2 = 0); the mean/variance formulas in .compute_sparse_pearson()
    # handle n = 0 correctly (the gene is then just a constant at `fill`).
    all_genes <- unname(row_keycol)
    missing_genes <- setdiff(all_genes, result$gene)
    if (length(missing_genes) > 0L) {
        result <- rbind(result, data.frame(gene = missing_genes, n = 0,
                                           sum_val = 0, S2 = 0))
    }

    result
}

#' @importFrom dplyr collect group_by summarize sql
#' @importFrom DuckDBDataFrame tblconn
.compute_pair_products <- function(table, row_key, col_key, datacol,
                                    pairings, row_keycol)
{
    # Use select = TRUE to get the transformed datacol values
    conn <- tblconn(table)

    # Build self-join query for the "both genes explicitly stored at this
    # cell" population (the only one of the 4 sparse sub-populations used
    # by .compute_sparse_pearson() that requires an actual join):
    #   SELECT a.gene AS gene1, b.gene AS gene2,
    #          COUNT(*) AS n11,
    #          SUM(a.value) AS sumA11, SUM(b.value) AS sumB11,
    #          COVAR_SAMP(a.value, b.value) AS cov11
    #   FROM data a JOIN data b ON a.cell = b.cell AND a.gene < b.gene
    #   GROUP BY a.gene, b.gene
    # COVAR_SAMP() is DuckDB's own numerically stable covariance aggregate
    # (never a hand-rolled SUM(a.value * b.value) raw cross-product).

    # Get underlying connection
    db_conn <- conn$src$con

    # Build SQL directly for self-join (dbplyr doesn't handle self-joins well)
    # Get the SQL for the base table
    base_sql <- dbplyr::sql_render(conn)
    datacol_name <- names(table@datacols)[1L]

    if (is.null(pairings)) {
        # All pairs
        query <- sprintf("
            SELECT a.\"%s\" as gene1, b.\"%s\" as gene2,
                   COUNT(*) as n11,
                   SUM(a.\"%s\") as sumA11,
                   SUM(b.\"%s\") as sumB11,
                   COVAR_SAMP(a.\"%s\", b.\"%s\") as cov11
            FROM (%s) a
            JOIN (%s) b ON a.\"%s\" = b.\"%s\" AND a.\"%s\" < b.\"%s\"
            GROUP BY a.\"%s\", b.\"%s\"
        ", row_key, row_key,
           datacol_name, datacol_name, datacol_name, datacol_name,
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
                   COUNT(*) as n11,
                   SUM(a.\"%s\") as sumA11,
                   SUM(b.\"%s\") as sumB11,
                   COVAR_SAMP(a.\"%s\", b.\"%s\") as cov11
            FROM (%s) a
            JOIN (%s) b ON a.\"%s\" = b.\"%s\"
            JOIN pairs p ON (a.\"%s\" = p.g1 AND b.\"%s\" = p.g2)
                         OR (a.\"%s\" = p.g2 AND b.\"%s\" = p.g1)
            WHERE a.\"%s\" < b.\"%s\"
            GROUP BY a.\"%s\", b.\"%s\"
        ", pairs_values,
           row_key, row_key,
           datacol_name, datacol_name, datacol_name, datacol_name,
           base_sql, base_sql,
           col_key, col_key,
           row_key, row_key, row_key, row_key,
           row_key, row_key,
           row_key, row_key)
    }

    result <- DBI::dbGetQuery(db_conn, query)
    result <- as.data.frame(result)

    # GOTCHA (same as elsewhere in this file): BIGINT COUNT(*) comes back
    # as bit64::integer64 from a raw DBI query too. Cast before arithmetic.
    result$n11 <- as.double(result$n11)

    result
}

#' @importFrom stats pt p.adjust
#
# .compute_sparse_pearson(): numerically-stable sparse Pearson correlation.
#
# HISTORY / WHY THE OLD FORMULA WAS UNSAFE:
# The original implementation computed, per gene, a raw-moment pair
# (sum(x), sum(x^2)) and, per gene pair, a raw cross-product sum
# sum(x * y) over the explicitly-stored cells, then recovered
# sum-of-squared-deviations and sum-of-cross-product-deviations via
#   SS  = sum(x^2) - sum(x)^2 / k
#   num = sum(x*y) - sum(x) * sum(y) / k
# This is the textbook "naive"/two-pass-avoiding formula for variance and
# covariance, and it is well known to be numerically unstable: when the
# data sit far from zero (large counts, pseudobulk sums, intensities, or
# just a nonzero `fill`), sum(x^2) and sum(x)^2/k are each huge and nearly
# equal, so their difference (the true, much smaller, sum of squared
# deviations) is computed as the difference of two large floating-point
# numbers -- catastrophic cancellation. In production this measurably
# produced |rho| > 1 and even NaN (from SS going NEGATIVE, which is
# mathematically impossible for a true sum of squared deviations) on
# realistic large-magnitude input. At small/log-normalized scale the old
# formula was fine, which is exactly what made this dangerous: it shipped
# silently wrong numbers only outside the regime the existing tests
# happened to exercise.
#
# THE STABLE REPLACEMENT:
# For a gene pair (A, B), partition the k cells into 4 disjoint
# populations by which of A/B has an EXPLICITLY STORED (non-fill) value
# at that cell:
#   11 = both stored      10 = A stored, B implicit (fill)
#   01 = A implicit, B stored      00 = both implicit (fill)
# The classic ANOVA-style decomposition of a cross-product sum into
# within-group + between-group pieces (the covariance generalization of
# the Chan et al. parallel-variance-combination identity; see
# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance) gives
#
#   Sum[(x - mu1)(y - mu2)] over k  =  Sum_g C_g
#                      + Sum_g n_g * (xbar_g - mu1) * (ybar_g - mu2)
#
# where C_g = Sum[(x - xbar_g)(y - ybar_g)] within group g, and mu1/mu2
# are the TRUE overall means over all k cells (stored + implicit), derived
# per-gene using the same stable Chan-et-al. fill-correction already used
# for modelGeneVar()/modelGeneVarByPoisson() (.fused_group_variance_stats()
# above): split into the n explicitly-stored cells (whose variance comes
# straight from DuckDB's native VAR_SAMP) and the (k - n) implicit cells
# pinned at the known constant `fill`.
#
# Three of the four C_g are EXACTLY ZERO: group 10 has gene B pinned at
# the constant `fill` for every cell in it, so gene B's within-group
# deviation is identically 0 there (and symmetrically for group 01's gene
# A); group 00 is pinned at `fill` for BOTH genes. Only C_11 -- the
# explicitly-stored/explicitly-stored population -- can be nonzero, and it
# comes straight from DuckDB's native COVAR_SAMP() over a real inner join
# (.compute_pair_products() above). No raw SUM(x * y) is ever computed,
# and no term below is ever a difference of two large near-equal numbers:
# every cross term is a count times a difference of two already-stable
# means, exactly the same algebraic shape already validated for the
# variance fix.
.compute_sparse_pearson <- function(row_stats, pair_stats,
                                     row_keycol, ncells, pairings, fill = 0)
{
    k <- as.numeric(ncells)

    # ---- Per-gene true mean/SS over ALL k cells (stored + implicit) ----
    rs <- row_stats
    rs$mu <- (rs$sum_val + fill * (k - rs$n)) / k
    ybar <- ifelse(rs$n > 0, rs$sum_val / rs$n, fill)
    rs$SS <- (rs$n - 1) * rs$S2 + rs$n * (ybar - rs$mu)^2 +
        (k - rs$n) * (fill - rs$mu)^2
    rs <- rs[, c("gene", "n", "sum_val", "mu", "SS")]

    # ---- Full pair skeleton, independent of which pairs the join query
    # ---- happened to return a row for (see below) ----
    if (is.null(pairings)) {
        genes <- unname(row_keycol)
        all_pairs <- expand.grid(gene1 = genes, gene2 = genes)
        all_pairs <- all_pairs[all_pairs$gene1 < all_pairs$gene2, , drop = FALSE]
    } else {
        if (is.character(pairings)) {
            pairings <- matrix(match(pairings, names(row_keycol)), ncol = 2)
        }
        keyed <- matrix(row_keycol[pairings], ncol = 2)
        all_pairs <- data.frame(gene1 = pmin(keyed[, 1], keyed[, 2]),
                                gene2 = pmax(keyed[, 1], keyed[, 2]))
    }

    # A gene pair with NO overlapping explicitly-stored cells at all (e.g.
    # fully disjoint supports) produces no row in the INNER JOIN behind
    # pair_stats. left-join onto the full pair skeleton and treat those as
    # n11 = 0 -- the group-11 population is simply empty for that pair,
    # which every formula below already handles without special-casing.
    pair_stats <- merge(all_pairs, pair_stats, by = c("gene1", "gene2"),
                        all.x = TRUE)
    pair_stats$n11[is.na(pair_stats$n11)] <- 0
    pair_stats$sumA11[is.na(pair_stats$sumA11)] <- 0
    pair_stats$sumB11[is.na(pair_stats$sumB11)] <- 0
    # cov11 is left NA for n11 < 2 (COVAR_SAMP is undefined there); every
    # use of cov11 below is guarded by n11 >= 2 so the NA never surfaces.

    pair_stats <- merge(pair_stats,
                        setNames(rs, c("gene1", "n1", "sum1", "mu1", "SS1")),
                        by = "gene1", all.x = TRUE)
    pair_stats <- merge(pair_stats,
                        setNames(rs, c("gene2", "n2", "sum2", "mu2", "SS2")),
                        by = "gene2", all.x = TRUE)

    # ---- Fold in the 10/01/00 populations analytically (see header) ----
    n11 <- pair_stats$n11
    n10 <- pair_stats$n1 - n11
    n01 <- pair_stats$n2 - n11
    n00 <- k - pair_stats$n1 - pair_stats$n2 + n11

    sumA10 <- pair_stats$sum1 - pair_stats$sumA11   # ybar10 == fill exactly
    sumB01 <- pair_stats$sum2 - pair_stats$sumB11   # xbar01 == fill exactly

    mu1 <- pair_stats$mu1
    mu2 <- pair_stats$mu2

    # C_11: the only possibly-nonzero within-group term. (n11 - 1) <= 0
    # whenever cov11 is undefined (n11 < 2), so guarding on n11 >= 2 keeps
    # that NA from ever being multiplied into the result.
    C11 <- ifelse(n11 >= 2, (n11 - 1) * pair_stats$cov11, 0)

    # Cross terms for groups 10/01/00, written as (group sum - n_g * mu)
    # rather than n_g * (group mean - mu) so no division by a possibly-zero
    # n_g ever appears (n_g * xbar_g IS the already-known group sum).
    term10 <- (sumA10 - n10 * mu1) * (fill - mu2)
    term01 <- (sumB01 - n01 * mu2) * (fill - mu1)
    term00 <- n00 * (fill - mu1) * (fill - mu2)

    # Group 11's cross term does need the actual group means; guarded on
    # n11 > 0 so 0/0 (when n11 == 0) never poisons the ifelse() result.
    xbar11 <- pair_stats$sumA11 / n11
    ybar11 <- pair_stats$sumB11 / n11
    term11 <- ifelse(n11 > 0, n11 * (xbar11 - mu1) * (ybar11 - mu2), 0)

    numerator <- C11 + term10 + term01 + term00 + term11
    denominator <- sqrt(pair_stats$SS1 * pair_stats$SS2)

    # Handle edge cases (zero variance)
    pair_stats$rho <- ifelse(denominator > 0, numerator / denominator, NA_real_)

    # Significance of the Pearson correlation via the standard t-test:
    # t = rho * sqrt((k-2) / (1-rho^2)), df = k-2. |rho| == 1 (a perfect,
    # degenerate correlation) makes the t-statistic formula divide by zero;
    # such a correlation is definitionally significant for k > 2 cells.
    df <- k - 2
    tstat <- pair_stats$rho * sqrt(df / (1 - pair_stats$rho^2))
    pair_stats$p.value <- ifelse(abs(pair_stats$rho) >= 1, 0, 2 * pt(-abs(tstat), df = df))
    pair_stats$p.value[df <= 0] <- NA_real_
    pair_stats$FDR <- p.adjust(pair_stats$p.value, method = "BH")

    # Return ordered by gene1, gene2
    pair_stats <- pair_stats[order(pair_stats$gene1, pair_stats$gene2),
                             c("gene1", "gene2", "rho", "p.value", "FDR"), drop = FALSE]
    rownames(pair_stats) <- NULL

    pair_stats
}
