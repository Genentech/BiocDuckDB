#' Options for DuckDB-backed matrices
#'
#' @description
#' Options for initializing DuckDB-backed matrices in
#' \code{\link{initializeCpp}}.
#'
#' @param option String specifying the name of the option.
#' @param value Value of the option.
#'
#' @details
#' The following options are supported:
#' \itemize{
#'   \item \code{realize}, a logical scalar specifying whether to materialize
#'     a \linkS4class{DuckDBArraySeed}'s underlying COO table into an
#'     in-memory sparse matrix with \code{\link{loadIntoMemory}} on every
#'     \code{\link{initializeCpp}} call (unlike \code{beachmat.hdf5}/
#'     \code{beachmat.tiledb}, this is deliberately \strong{not} cached across
#'     calls, see \code{\link{initializeCpp}}'s details for why). Unlike
#'     \code{beachmat.hdf5}/\code{beachmat.tiledb} (which default this to
#'     \code{FALSE} because they have a genuine on-demand streaming backend),
#'     this defaults to \code{TRUE}: BiocDuckDB has no native out-of-core
#'     \pkg{tatami} backend, so setting this to \code{FALSE} falls back to
#'     \pkg{beachmat}'s generic "unknown matrix" block-processing path rather
#'     than a faster streaming one.
#'   \item \code{memory_limit}, a numeric scalar specifying the maximum
#'     estimated size, in bytes, of the sparse matrix materialized by
#'     \code{\link{loadIntoMemory}}. Above this limit, materialization is
#'     refused so a genuinely out-of-core matrix fails loudly (and falls back
#'     to the generic path, when reached via \code{\link{initializeCpp}})
#'     instead of risking an OOM.
#' }
#'
#' @return If \code{value} is missing, the current setting of \code{option} is
#' returned. If \code{value} is supplied, it is used to set the option, and
#' the previous value of the option is invisibly returned.
#'
#' @author Patrick Aboyoun
#'
#' @seealso
#' \code{\link{initializeCpp}}, \code{\link{loadIntoMemory}}
#'
#' @examples
#' initializeOptions("memory_limit")
#' old <- initializeOptions("memory_limit", 4 * 1024^3)
#' initializeOptions("memory_limit")
#' initializeOptions("memory_limit", old)
#'
#' @export
#' @name initializeOptions
initializeOptions <- function(option, value) {
    old <- get(option, envir = .duckdb_beachmat_options, inherits = FALSE)
    if (missing(value)) {
        return(old)
    }
    assign(option, value, envir = .duckdb_beachmat_options)
    invisible(old)
}

.duckdb_beachmat_options <- new.env()
.duckdb_beachmat_options$realize <- TRUE
.duckdb_beachmat_options$memory_limit <- 2 * 1024^3
