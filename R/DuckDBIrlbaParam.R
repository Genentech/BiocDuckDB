#' A fast SVD path for DuckDB-backed matrices
#'
#' @description
#' \code{DuckDBIrlbaParam} is a \code{\link[BiocSingular]{BiocSingularParam}}
#' for \code{\link[BiocSingular]{runSVD}} that, for a \linkS4class{DuckDBMatrix},
#' materializes the matrix once into an in-memory sparse matrix and drives
#' \pkg{irlba} directly on it, instead of letting \pkg{irlba}'s Lanczos
#' iteration call \code{\%*\%} on the lazy \linkS4class{DuckDBMatrix} once per
#' solver step.
#'
#' @details
#' \code{\link[DuckDBArray]{DuckDBMatrix}} already implements \code{\%*\%}/
#' \code{crossprod} as correct SQL-pushdown methods (see \pkg{DuckDBArray}), so
#' \code{BiocSingular::runSVD(x, BSPARAM = BiocSingular::IrlbaParam())} already
#' computes a correct PCA/SVD directly on a \code{DuckDBMatrix}. But each of
#' those \code{\%*\%} calls pays real SQL query-construction cost, and
#' \pkg{irlba}'s Lanczos loop calls it once per solver iteration, hundreds of
#' times per solve, so the ordinary \code{IrlbaParam} path is markedly slower
#' on a \code{DuckDBMatrix} than the same computation in memory (measured: on
#' a real 12,500-cell x 200-HVG benchmark, ~16s vs ~0.3s in-memory).
#'
#' \code{DuckDBIrlbaParam} closes that gap for matrices that fit in memory, by
#' reusing the same \code{\link{loadIntoMemory}} machinery
#' (\code{.duckdb_seed_to_sparse_matrix()}, size-gated by
#' \code{"memory_limit"}) to materialize the matrix exactly once, then calling
#' \pkg{BiocSingular}'s own exported \code{irlba}-driving function
#' (\code{\link[BiocSingular]{runIrlbaSVD}}) directly on the materialized
#' matrix. Reusing that function, rather than reimplementing its nu/nv
#' trimming, edge cases, and output packaging, avoids a larger and riskier
#' undertaking than depending on one stable, public function.
#'
#' Whenever the fast path is not applicable (\code{x} is not a
#' \linkS4class{DuckDBMatrix}, its seed is not zero-filled, or the estimated
#' materialized size is over \code{"memory_limit"}), this falls back to
#' ordinary \code{IrlbaParam} behavior (the lazy SQL-pushdown \code{\%*\%}
#' path), so it can never turn previously-\emph{working} (if slower) code into
#' a failure. It can still fail the same way the ordinary path already does,
#' e.g. a non-zero-filled seed is not supported by \code{\%*\%} either, that
#' case falls back to, and fails identically to, the ordinary path rather than
#' silently producing an incorrect result.
#'
#' @param memory_limit Numeric scalar, see the \code{memory_limit} option in
#'   \code{\link{initializeOptions}}. Materialization is refused above this
#'   estimated size (bytes), falling back to the ordinary lazy path.
#' @param deferred,fold,extra.work,... Passed to
#'   \code{\link[BiocSingular]{IrlbaParam}}.
#' @param x A matrix-like object, typically a \linkS4class{DuckDBMatrix} (the
#'   fast path only applies to one; any other input falls back to ordinary
#'   \code{IrlbaParam} behavior).
#' @param k,nu,nv,center,scale,BPPARAM,BSPARAM See
#'   \code{\link[BiocSingular]{runSVD}}.
#'
#' @return A \code{DuckDBIrlbaParam} object, for use as the \code{BSPARAM}
#' argument to \code{\link[BiocSingular]{runSVD}} (and so
#' \code{scran::fixedPCA()}, \code{scater::runPCA()}/\code{calculatePCA()}).
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[BiocSingular]{IrlbaParam}} for the ordinary, lazy-pushdown behavior
#'   \item \code{\link{loadIntoMemory}} for the materialization helper this reuses
#'   \item \code{\link{initializeOptions}} for the shared \code{"memory_limit"} default
#' }
#'
#' @examples
#' df <- data.frame(row = c("r1", "r2", "r1"), col = c("c1", "c1", "c2"),
#'                  value = c(1, 2, 3))
#' tf <- tempfile(fileext = ".parquet")
#' on.exit(unlink(tf))
#' arrow::write_parquet(df, tf)
#' pqmat <- DuckDBArray::DuckDBMatrix(tf, datacol = "value",
#'                                    keycols = list(row = c("r1", "r2"), col = c("c1", "c2")))
#' BiocSingular::runSVD(pqmat, k = 1, BSPARAM = DuckDBIrlbaParam())
#'
#' @aliases DuckDBIrlbaParam-class
#' @aliases runSVD,DuckDBIrlbaParam-method
#'
#' @include initializeCpp.R
#' @include initializeOptions.R
#'
#' @name DuckDBIrlbaParam
NULL

#' @export
#' @importClassesFrom BiocSingular IrlbaParam
setClass("DuckDBIrlbaParam", contains = "IrlbaParam",
         slots = c(memory_limit = "numeric"))

#' @export
#' @rdname DuckDBIrlbaParam
DuckDBIrlbaParam <-
function(memory_limit = initializeOptions("memory_limit"),
         deferred = FALSE, fold = Inf, extra.work = 7, ...)
{
    new("DuckDBIrlbaParam", deferred = as.logical(deferred),
        fold = as.numeric(fold), extra.work = as.integer(extra.work),
        args = list(...), memory_limit = as.numeric(memory_limit))
}

#' @export
#' @rdname DuckDBIrlbaParam
#' @importFrom BiocSingular ExactParam runIrlbaSVD runSVD
#' @importFrom BiocParallel SerialParam
#' @importClassesFrom DuckDBArray DuckDBMatrix
setMethod("runSVD", "DuckDBIrlbaParam", function(
    x, k, nu = k, nv = k, center = FALSE, scale = FALSE,
    BPPARAM = SerialParam(), ..., BSPARAM)
{
    result <- tryCatch({
        if (!is(x, "DuckDBMatrix")) {
            stop("not a DuckDBMatrix")
        }
        sp <- .duckdb_seed_to_sparse_matrix(x@seed, BSPARAM@memory_limit)
        do.call(runIrlbaSVD,
            c(list(sp, k = k, nu = nu, nv = nv, center = center, scale = scale,
                   deferred = BSPARAM@deferred, extra.work = BSPARAM@extra.work,
                   fold = BSPARAM@fold, BPPARAM = BPPARAM), BSPARAM@args))
    }, error = function(e) e)
    if (!inherits(result, "error")) {
        return(result)
    }
    callNextMethod()
})
