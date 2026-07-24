#' Stream a large table to a single multi-part Parquet resource
#'
#' @description
#' Write a table that does not fit in memory to one flat, multi-part Parquet
#' resource by pulling it block-by-block from a producer callback. This owns the
#' streaming bookkeeping that every large producer would otherwise hand-roll on
#' top of \code{\link{writeParquet}}: the running row \code{offset}, the
#' \code{part} index and zero-padded \code{part_digits}, the \code{append} flag
#' (part 0 creates the directory, later parts append), the positional
#' \code{dimtbl} slice per block, the \code{index_max} threading that keeps the
#' \code{__index__} column one consistent (possibly >2^31) integer type across
#' parts, and the post-write integrity checks. The single descriptor it returns
#' is ready to hand to \code{\link{writeDatapackage}}.
#'
#' @param blocks A function of one argument, the 1-based part index \code{i}. It
#'   must return the \code{i}-th block as a \code{data.frame} (or
#'   \code{DataFrame} / anything \code{as.data.frame}-able), or \code{NULL} when
#'   the stream is exhausted. Blocks are consumed in ascending index order; how
#'   a block maps to source rows (a row range, a sample id, a query page) is the
#'   producer's business. Empty (zero-row) blocks are skipped without consuming
#'   a \code{part}.
#' @param path Directory for the resource's Parquet parts.
#' @param dimension Biological axis of the resource: one of \code{"unbound"},
#'   \code{"sample"}, \code{"feature"}, \code{"crossed"}.
#' @param indexcol,keycol Passed to \code{\link{writeParquet}}. \code{indexcol}
#'   is the streamed row index (default \code{"__index__"}); \code{keycol} is
#'   the row-name column (\code{"__name__"}), or \code{NULL} for a row-number
#'   key.
#' @param layout Physical layout string (default \code{"data_frame"}); also
#'   accepts other flat layouts \code{\link{writeParquet}} supports for a
#'   \code{data.frame} (e.g. \code{"spatial_points"}).
#' @param refs,cluster_by Passed to \code{\link{writeParquet}} (foreign-key
#'   references; on-write clustering key).
#' @param dimtbl An optional full, index-ordered dimension-table
#'   \code{DataFrameList} for the resource. It is sliced positionally per block
#'   (\code{dimtbl[offset + seq_len(nrow(block)), ]}), so it must have exactly
#'   one row per streamed row in the same order. A drift is fatal (see Details).
#' @param index_max Optional upper bound on the row index (the resource's total
#'   row count, or \code{Inf} when unknown ahead of streaming). Passing it types
#'   \code{__index__} wide enough on part 0 for the whole resource, which is
#'   \strong{required} for a > 2^31-row resource and removes any dependence on
#'   block size (see Details). \code{NULL} keeps per-block index narrowing.
#' @param expected_rows Optional expected total row count. When given and the
#'   streamed total differs, a \emph{warning} is emitted (an unbound junction
#'   table can legitimately differ; a bound one usually indicates an off-by-one
#'   or a dropped block).
#' @param part_digits Zero-padding width for the part index in filenames. The
#'   reader globs \code{*.parquet}, so this is cosmetic; defaults to \code{6L}.
#' @param name Resource name recorded in the descriptor (default
#'   \code{basename(path)}).
#' @param ... Further arguments forwarded to \code{\link{writeParquet}}.
#'
#' @details
#' \strong{Index typing / the narrowing floor.} \code{\link{writeParquet}} types
#' part 0's \code{__index__} column by the index range it sees. Without
#' \code{index_max}, part 0 is typed to its own \code{[0, nrow]} range, so a
#' small part 0 can pick a narrow integer type (e.g. \code{uint16}) that a later
#' part's larger index overflows on append. Passing \code{index_max} (the total
#' row count, or \code{Inf}) types the index for the whole resource up front and
#' removes this hazard entirely. When \code{index_max} is not given and part 0
#' has fewer than 65536 rows in a multi-part stream, a warning is emitted.
#'
#' \strong{Partition alignment.} The per-block \code{dimtbl} slice assumes the
#' stream is dense and contiguous with exactly \code{nrow(dimtbl)} rows in index
#' order. If the streamed total does not equal \code{nrow(dimtbl)} the partition
#' column would be misaligned (silently corrupting directory pruning), so this
#' is a fatal error rather than a warning.
#'
#' @return
#' Invisibly, the single Frictionless resource descriptor (a list) captured from
#' part 0, with an added \code{n_rows} count, ready for
#' \code{\link{writeDatapackage}}. \code{NULL} if the stream yielded no rows.
#'
#' @seealso \code{\link{writeParquet}} for the single-block writer;
#'   \code{\link{writeDatapackage}} to assemble the resources into a
#'   \code{datapackage.json}.
#'
#' @examples
#' tf <- tempfile()
#' on.exit(unlink(tf, recursive = TRUE))
#' chunks <- list(data.frame(v = 1:3), data.frame(v = 4:5))
#' res <- writeStreamingResource(
#'     function(i) if (i <= length(chunks)) chunks[[i]] else NULL,
#'     path = file.path(tf, "samples"), dimension = "unbound",
#'     keycol = NULL, index_max = 5, expected_rows = 5)
#' res[[1]][["n_rows"]]                  # 5
#' list.files(file.path(tf, "samples"))  # part-000000.parquet, part-000001.parquet
#'
#' @author Patrick Aboyoun
#'
#' @include writeParquet.R
#'
#' @export
writeStreamingResource <-
function(blocks, path, dimension,
         indexcol = "__index__", keycol = "__name__", layout = "data_frame",
         refs = NULL, dimtbl = NULL, cluster_by = NULL,
         index_max = NULL, expected_rows = NULL,
         part_digits = 6L, name = basename(path), ...)
{
    if (!is.function(blocks)) {
        stop("'blocks' must be a function of one argument (the 1-based part ",
             "index) returning a data.frame block, or NULL when exhausted")
    }

    resource <- NULL
    offset <- 0            # numeric (not 0L): a > 2^31-row resource overflows
    written <- 0L          # count of NON-EMPTY parts actually written
    part0_rows <- NA_real_ # rows in part 0, for the narrowing-floor warning
    i <- 0L
    repeat {
        i <- i + 1L
        block <- blocks(i)
        if (is.null(block)) {
            break
        }
        if (!is.data.frame(block)) {
            block <- as.data.frame(block, optional = TRUE)
        }
        if (nrow(block) == 0L) {
            next   # skip empty blocks so part 0 is the first block with rows
        }

        # Narrowing-floor warning: without index_max, part 0's __index__ is
        # typed to its own [0, nrow] range; a later part's larger index then
        # overflows that narrowed type on append. Warn when a small part 0 is
        # followed by more parts and no index_max was supplied.
        if (written == 1L && is.null(index_max) &&
            !is.na(part0_rows) && part0_rows < 65536) {
            warning(sprintf(paste0(
                "%s: part 0 has %.0f (< 65536) rows and 'index_max' was not set; ",
                "appending further parts may overflow the narrowed '%s' integer ",
                "type. Pass index_max = <total rows> for a multi-part stream."),
                basename(path), part0_rows, indexcol), call. = FALSE)
        }

        # Positional dimtbl slice for this block (index-ordered, offset-aligned)
        dimtbl_block <- NULL
        if (!is.null(dimtbl)) {
            dimtbl_block <- dimtbl[offset + seq_len(nrow(block)), , drop = FALSE]
        }

        res <- writeParquet(block, path = path, indexcol = indexcol,
                            keycol = keycol, dimension = dimension,
                            layout = layout, refs = refs, dimtbl = dimtbl_block,
                            cluster_by = cluster_by, offset = offset,
                            part = written, part_digits = part_digits,
                            append = written > 0L, index_max = index_max,
                            name = name, ...)
        if (written == 0L) {
            resource <- res      # descriptor is returned only from part 0
            part0_rows <- nrow(block)
        }
        written <- written + 1L
        offset <- offset + nrow(block)
    }

    # Partition-alignment guard: a per-block dimtbl slice assumes a dense,
    # contiguous stream of exactly nrow(dimtbl) rows. A drift writes NA group
    # labels, silently corrupting directory pruning -- fail loudly.
    if (!is.null(dimtbl) && offset != nrow(dimtbl)) {
        stop(sprintf(paste0("streamed %.0f rows but dimtbl has %d -- the ",
                            "partition column would be misaligned"),
                     offset, nrow(dimtbl)))
    }

    # Coverage check: a declared expected row count that the stream did not
    # deliver usually means a dropped block or an off-by-one (warning: an
    # unbound junction table can legitimately differ).
    if (!is.null(expected_rows) && offset != expected_rows) {
        warning(sprintf(paste0("%s streamed %.0f rows but %d were expected ",
                               "(%+.0f)."), basename(path), offset,
                        expected_rows, offset - expected_rows), call. = FALSE)
    }

    # Stamp the true row count onto the descriptor for downstream checks.
    if (!is.null(resource) && length(resource)) {
        resource[[length(resource)]][["n_rows"]] <- offset
    }
    invisible(resource)
}
