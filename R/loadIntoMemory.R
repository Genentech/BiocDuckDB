#' Load a DuckDB-backed matrix into memory
#'
#' @description
#' Load a \linkS4class{DuckDBMatrix}'s (or its \linkS4class{DuckDBArraySeed}'s)
#' underlying COO table into memory as an external pointer to a
#' \pkg{tatami}-compatible representation, via a single bulk SQL pull. This is
#' what \code{\link{initializeCpp}} calls internally when
#' \code{duckdb.realize = TRUE} (the default); call it directly when you want
#' the materialized pointer without going through the \code{initializeCpp}
#' generic. Each call re-materializes; the result is deliberately not cached
#' (see \code{\link{initializeCpp}}'s details for why).
#'
#' @param x A \linkS4class{DuckDBMatrix} or \linkS4class{DuckDBArraySeed}
#'   object.
#' @param memory_limit Numeric scalar, see the \code{memory_limit} option in
#'   \code{\link{initializeOptions}}. Materialization is refused above this
#'   estimated size (bytes) rather than risking an OOM.
#'
#' @return An external pointer that can be used in \pkg{tatami}-based functions.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \code{\link{initializeCpp}}, \code{\link{initializeOptions}}
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
#' ptr <- loadIntoMemory(pqmat)
#'
#' @export
#' @importClassesFrom DuckDBArray DuckDBArraySeed DuckDBMatrix
#' @importFrom beachmat initializeCpp
#' @include initializeCpp.R
loadIntoMemory <-
function(x, memory_limit = initializeOptions("memory_limit"))
{
    if (is(x, "DuckDBMatrix")) {
        x <- x@seed
    }
    if (!is(x, "DuckDBArraySeed")) {
        stop("unsupported seed type '", class(x)[1L], "'")
    }
    initializeCpp(.duckdb_seed_to_sparse_matrix(x, memory_limit))
}
