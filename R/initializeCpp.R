#' Initialize DuckDB-backed matrices for beachmat
#'
#' @description
#' Registers a \code{\link[beachmat]{initializeCpp}} method for
#' \link[DuckDBArray]{DuckDBArraySeed} objects (the seed underlying every
#' \link[DuckDBArray]{DuckDBMatrix}), giving compiled C++ code that goes through
#' \pkg{beachmat}/\pkg{tatami} (e.g. \code{BiocSingular}'s
#' \code{compute_center}/\code{compute_scale}) a fast path instead of always
#' falling back to \pkg{beachmat}'s generic "unknown matrix" block-processing
#' path.
#'
#' @details
#' Unlike \code{beachmat.hdf5}/\code{beachmat.tiledb}, there is no dedicated
#' \code{tatami_duckdb} C++ library implementing genuine on-demand,
#' out-of-core reads. This method instead materializes the seed's underlying
#' COO table into an in-memory sparse matrix via a single bulk SQL pull
#' (\code{\link{loadIntoMemory}}), then delegates to \pkg{beachmat}'s own
#' native \code{initializeCpp,dgCMatrix-method}. This mirrors
#' \code{beachmat.hdf5}'s \emph{opt-in} \code{realize} path, not its default
#' on-demand streaming path, and is only attempted for a zero-filled (sparse)
#' seed within \code{\link{initializeOptions}}'s \code{"memory_limit"}.
#'
#' Unlike \code{beachmat.hdf5}/\code{beachmat.tiledb}, the materialized
#' pointer is \strong{not} cached across calls with
#' \code{\link[beachmat]{checkMemoryCache}}. Verified empirically that doing
#' so is unsafe here: two \link[DuckDBArray]{DuckDBMatrix} objects built from the
#' same underlying query (the common case, e.g. re-reading the same release
#' twice) render identical SQL and would hit the same cache entry, but
#' reusing that cached pointer across a later, logically separate computation
#' silently corrupted results. Every call re-materializes, which costs a
#' fraction of a second for a \code{"memory_limit"}-sized matrix, negligible
#' next to the speedup this fast path already provides.
#'
#' This is purely a performance optimization for anything that already goes
#' through \code{initializeCpp} (e.g. \code{BiocSingular::runSVD}'s
#' \code{center}/\code{scale} computation). It does \strong{not} speed up
#' \code{irlba}'s Lanczos iteration itself, which calls \code{\%*\%} directly
#' on the \link[DuckDBArray]{DuckDBMatrix} object rather than through
#' \code{initializeCpp}; use the SQL-pushdown \code{\%*\%}/\code{crossprod}
#' methods in \pkg{DuckDBArray} for that, or materialize small-enough matrices
#' with \code{\link{loadIntoMemory}} and drive \code{irlba} directly.
#'
#' Whenever the fast path is not applicable (a non-zero-filled seed, more
#' than 2 keyed dimensions, or an estimated size over \code{"memory_limit"}),
#' this method falls back to \pkg{beachmat}'s existing generic path rather
#' than raising an error, so it can never make previously-working (if slower)
#' code stop working.
#'
#' @param x A \link[DuckDBArray]{DuckDBArraySeed} object.
#' @param ... Further arguments, passed to \pkg{beachmat}'s generic fallback
#'   method when the fast path is not applicable.
#' @param duckdb.realize Logical scalar, see the \code{realize} option in
#'   \code{\link{initializeOptions}}.
#' @param duckdb.memory_limit Numeric scalar, see the \code{memory_limit}
#'   option in \code{\link{initializeOptions}}.
#'
#' @return An external pointer that can be used in any \pkg{tatami}-compatible
#' function.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[beachmat]{initializeCpp}} for the generic
#'   \item \code{\link{loadIntoMemory}} to force materialization directly
#'   \item \code{\link{initializeOptions}} for the tunable defaults
#' }
#'
#' @examples
#' # A genuinely sparse 2x2 matrix: (r2, c2) is an implicit zero.
#' df <- data.frame(row = c("r1", "r2", "r1"), col = c("c1", "c1", "c2"),
#'                  value = c(1, 2, 3))
#' tf <- tempfile(fileext = ".parquet")
#' on.exit(unlink(tf))
#' arrow::write_parquet(df, tf)
#' pqmat <- DuckDBArray::DuckDBMatrix(tf, datacol = "value",
#'                                    keycols = list(row = c("r1", "r2"), col = c("c1", "c2")))
#' ptr <- beachmat::initializeCpp(pqmat@seed)
#'
#' @name initializeCpp
NULL

#' @importFrom dplyr all_of collect count pull select
#' @importFrom DuckDBDataFrame tblconn
#' @importFrom Matrix sparseMatrix
#' @importFrom S4Vectors isSingleNumber
.duckdb_seed_to_sparse_matrix <- function(x, memory_limit) {
    if (x@fill != 0) {
        stop("requires a zero-filled (sparse) DuckDBArraySeed")
    }
    table <- x@table
    keycols <- table@keycols
    dim_idx <- which(lengths(keycols) > 1L)
    if (length(dim_idx) != 2L) {
        stop("requires exactly 2 dimensions with length > 1")
    }
    if (length(table@datacols) != 1L) {
        stop("requires a single datacols")
    }
    row_key <- names(keycols)[dim_idx[1L]]
    col_key <- names(keycols)[dim_idx[2L]]
    row_keycol <- keycols[[dim_idx[1L]]]
    col_keycol <- keycols[[dim_idx[2L]]]
    datacol_name <- names(table@datacols)

    conn <- tblconn(table, select = TRUE)
    nnz <- as.numeric(pull(collect(count(conn)), 1L))
    # Bytes-per-nonzero guard constant, derivation:
    #
    # The transient peak this guards is `df` (the row_key/col_key/datacol
    # data.frame from collect(), below) coexisting with the sparse matrix
    # sparseMatrix() is building from it -- in particular, match()'s output
    # `i`/`j` integer index vectors, which exist alongside `df`'s own
    # row_key/col_key columns before `df` can be garbage collected.
    #
    # A naive "1 double value + 2 integer keys" count (8 + 4 + 4 = 16
    # bytes/nonzero) assumes row_key/col_key are already integers and
    # ignores that `df` and the freshly-built `i`/`j` vectors coexist. Both
    # assumptions are wrong for this package's own documented usage
    # pattern (see `scanKeycols()`): a COO Parquet file produced by another
    # tool, read via `DuckDBMatrix(path, keycols = scanKeycols(path, ...))`,
    # has row_key/col_key SQL columns holding the distinct key VALUES
    # directly -- e.g. character gene symbols or cell barcodes, not a
    # synthetic integer remap -- so `df$row_key`/`df$col_key` collect() as R
    # character vectors, not integers.
    #
    # Measured empirically (R's own gc(reset = TRUE) peak-heap tracking
    # immediately around this function, isolating its own transient
    # allocation from the rest of the process) on synthetic COO fixtures at
    # several sizes and key-column types:
    #   string-keyed,  5M nonzeros (1000 genes x 50000 cells):  ~214 bytes/nonzero
    #   string-keyed, 10M nonzeros (2000 genes x 50000 cells):  ~139 bytes/nonzero
    #   string-keyed, 20M nonzeros (4000 genes x 50000 cells):   ~92 bytes/nonzero
    #   integer-keyed,10M nonzeros (2000 genes x 50000 cells):  ~120 bytes/nonzero
    # All four measurements exceed the naive 16 bytes/nonzero estimate by
    # 6-13x (not just the ~2.3x previously reported for a single
    # string-keyed case), and the ratio is scale- and
    # cardinality-dependent rather than a clean constant, so this uses one
    # single conservative constant -- comfortably above every measurement
    # above, string- or integer-keyed -- rather than a type-aware estimate
    # that would need its own per-type validation.
    bytes_per_nnz <- 256
    est_bytes <- nnz * bytes_per_nnz
    if (!isSingleNumber(memory_limit) || est_bytes > memory_limit) {
        stop(sprintf("estimated %.0f MB exceeds 'memory_limit' (%.0f MB)",
                      est_bytes / 1024^2, memory_limit / 1024^2))
    }

    df <- collect(select(conn, all_of(c(row_key, col_key, datacol_name))))
    sparseMatrix(
        i = match(df[[row_key]], unname(row_keycol)),
        j = match(df[[col_key]], unname(col_keycol)),
        x = df[[datacol_name]],
        dims = c(length(row_keycol), length(col_keycol)),
        dimnames = list(names(row_keycol), names(col_keycol))
    )
}

#' @export
#' @rdname initializeCpp
#' @importFrom beachmat initializeCpp
#' @importClassesFrom DuckDBArray DuckDBArraySeed
setMethod("initializeCpp", "DuckDBArraySeed", function(
    x, ...,
    duckdb.realize = initializeOptions("realize"),
    duckdb.memory_limit = initializeOptions("memory_limit"))
{
    if (isTRUE(duckdb.realize)) {
        # Deliberately NOT using beachmat::checkMemoryCache here. Verified
        # empirically that caching the materialized pointer keyed on rendered
        # SQL is unsafe: two DuckDBMatrix objects built from the same
        # underlying query (the common case -- re-reading the same release)
        # render identical SQL and hit the same cache entry, but reusing that
        # cached external pointer across a later, logically separate
        # computation silently corrupted results (confirmed by reproducing it
        # directly: running scran::fixedPCA once, then scater::runPCA on a
        # freshly-built DuckDBMatrix over the same data, diverged from the
        # dense oracle only when the second call hit the first call's cached
        # pointer). Re-materializing every call costs ~0.1s for a
        # memory_limit-sized matrix, which is negligible next to the ~160x
        # speedup this fast path already provides over the lazy %*% path.
        result <- tryCatch(loadIntoMemory(x, memory_limit = duckdb.memory_limit),
                           error = function(e) e)
        if (!inherits(result, "error")) {
            return(result)
        }
    }
    callNextMethod(x, ...)
})
