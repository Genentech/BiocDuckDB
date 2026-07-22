### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Datapackage envelope assembler
###
### The public, source-agnostic seam between "write resources" and "write the
### manifest". The experiment-level writeParquet() methods accumulate a list of
### Frictionless resource descriptors and then assemble the datapackage.json
### envelope; this function is that envelope step, factored out so that (a) the
### SummarizedExperiment and MultiAssayExperiment methods single-source it and
### (b) producers that build resources piecemeal -- streaming a dataset too large
### to hold in memory, or promoting from a foreign store -- can emit a conformant
### datapackage.json without hand-rolling it. It is the write half of the
### ingest contract; the read half is readParquet() (container level) and the
### DuckDBMatrix()/DuckDBArray()/DuckDBTable() constructors (component level),
### which attach existing Parquet in place.

#' Write a Frictionless datapackage.json envelope
#'
#' Assembles and writes the top-level \code{datapackage.json} manifest from a
#' list of already-written Frictionless resource descriptors. This is the
#' envelope step shared by the experiment-level \code{\link{writeParquet}}
#' methods, exposed so that producers who accumulate resources incrementally
#' (streaming large datasets, or promoting from another store) can emit a
#' conformant manifest without reconstructing an in-memory Bioconductor object.
#'
#' Each entry of \code{resources} is a Frictionless resource descriptor -- a list
#' with \code{name}, \code{path}, \code{dimension}, \code{layout}, \code{format},
#' \code{mediatype}, and \code{schema} -- exactly as returned by the primitive
#' \code{\link{writeParquet}} methods (array, \code{data.frame}, \code{DataFrame},
#' \code{SelfHits}). \code{NULL} entries are dropped, so the \code{NULL} returned
#' by append/streaming parts (see \code{\link{writeParquet}}) can be accumulated
#' and passed straight through. Descriptors are written verbatim otherwise; strip
#' any private, non-Frictionless keys before calling.
#'
#' @param model Character(1) package-level schema identifier that selects the
#'   \code{\link{readParquet}} reader used to reconstruct the container (e.g.
#'   \code{"summarized_experiment"}, \code{"single_cell_experiment"},
#'   \code{"multi_assay_experiment"}). See the storage-layout vignette for the
#'   documented \code{model} values.
#' @param resources A list of Frictionless resource descriptors (each a list),
#'   as returned/accumulated from \code{\link{writeParquet}}. \code{NULL} entries
#'   are removed.
#' @param path Character(1) directory to write \code{datapackage.json} into;
#'   created recursively if it does not exist.
#' @param main_exp_name Optional character(1) naming the main experiment (used by
#'   the \code{single_cell_experiment} reader). Omitted from the manifest when
#'   \code{NULL}.
#' @param annotations Optional list of non-relational metadata elements (as
#'   produced during metadata serialization). Omitted when \code{NULL}.
#'
#' @return Invisibly, the assembled package list that was written.
#'
#' @examples
#' # Assemble a manifest from a hand-built resource descriptor.
#' tf <- tempfile()
#' resources <- list(list(
#'     name = "features", path = "features",
#'     dimension = "feature", layout = "data_frame",
#'     format = "parquet",
#'     mediatype = "application/vnd.apache.parquet",
#'     schema = list(fields = list(list(name = "id", type = "integer")))))
#' writeDatapackage("summarized_experiment", resources, tf)
#' cat(readLines(file.path(tf, "datapackage.json")), sep = "\n")
#'
#' @seealso \code{\link{writeParquet}} for writing resources, and
#'   \code{\link{readParquet}} for reading a written package (the
#'   \code{DuckDBMatrix}/\code{DuckDBArray}/\code{DuckDBTable} constructors attach
#'   an existing coord-array in place).
#'
#' @author Patrick Aboyoun
#'
#' @importFrom jsonlite write_json
#' @importFrom S4Vectors isSingleString
#' @export
writeDatapackage <- function(model, resources, path,
                             main_exp_name = NULL, annotations = NULL)
{
    if (!isSingleString(model))
        stop("'model' must be a single non-NA string")
    if (!isSingleString(path))
        stop("'path' must be a single non-NA string")
    if (!is.list(resources))
        stop("'resources' must be a list of Frictionless resource descriptors")
    if (!is.null(main_exp_name) && !isSingleString(main_exp_name))
        stop("'main_exp_name' must be NULL or a single string")

    # Drop NULL descriptors (append/streaming parts return NULL).
    resources <- Filter(Negate(is.null), resources)

    package <- list(model = model, resources = resources)
    if (!is.null(main_exp_name))
        package[["main_exp_name"]] <- main_exp_name
    if (!is.null(annotations))
        package[["annotations"]] <- annotations

    # Declare the Frictionless profile as the leading key.
    package <- c(list("$schema" = .BIOCDUCKDB_PROFILE), package)

    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    write_json(package, path = file.path(path, "datapackage.json"),
               auto_unbox = TRUE, pretty = TRUE)

    invisible(package)
}
