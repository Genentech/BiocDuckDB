# shared streaming-output READ contract, asserted R-side.
#
# ratifies engine-local write *layout* under a *shared read contract*. These
# tests encode that contract's physical invariants as named assertions on
# BiocDuckDB-written streamed output.
#
# 1.  glob-readable: a resource is a directory of part-*.parquet read as a
#     whole.
# 2.  sorted by the key: per-row-group min/max zonemaps are non-overlapping.
# 3.  row-groups-per-file scale via a FIXED row-group size (491520), so a large
#     resource always has enough groups to saturate query threads.
# 4.  per-row-group min/max statistics are present (the pruning substrate).
# 5.  commit boundary: a completed stream is a clean parquet-only directory; an
#     interrupted one keeps an _INCOMPLETE marker and is not silently readable.

library(S4Vectors)

# --- a streamed multi-part resource reads back whole + in order ---
test_that("streamed multi-part resource is glob-readable and key-ordered", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    chunks <- list(
        data.frame(v = 1:3, g = c("x", "y", "z"), stringsAsFactors = FALSE),
        data.frame(v = 4:5, g = c("p", "q"), stringsAsFactors = FALSE))
    res <- writeStreamingResource(
        function(i) if (i <= length(chunks)) chunks[[i]] else NULL,
        path = file.path(tf, "tbl"), dimension = "unbound", keycol = NULL,
        index_max = 5, expected_rows = 5, name = "tbl")

    # genuinely multi-part (sequential part files), the whole set globbed.
    expect_length(list.files(file.path(tf, "tbl"), pattern = "\\.parquet$"), 2L)

    # commit boundary: a completed stream removes its in-progress marker, so a
    # finished resource is a clean directory of only parquet parts.
    expect_false(file.exists(file.path(tf, "tbl", "_INCOMPLETE")))

    tbl <- DuckDBDataFrame::DuckDBDataFrame(file.path(tf, "tbl"),
                                            keycol = "__index__")
    idx <- as.integer(rownames(as.data.frame(tbl)))
    # INV-2: the key is contiguous + non-decreasing across parts (zonemaps prune).
    expect_identical(sort(idx), 1:5)
    expect_false(is.unsorted(sort(idx)))
    expect_identical(res[[1L]][["n_rows"]], 5)
})

# --- an interrupted stream keeps the marker and is not silently readable ---
test_that("an incomplete resource (marker present) fails loudly, not partial", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    writeStreamingResource(
        function(i) if (i == 1L) data.frame(v = 1:3) else NULL,
        path = file.path(tf, "tbl"), dimension = "unbound", keycol = NULL,
        index_max = 3, name = "tbl")
    file.create(file.path(tf, "tbl", "_INCOMPLETE"))   # simulate a crashed write

    # The reader refuses the marked directory explicitly. This assertion used
    # to be a bare expect_error(), which passed for the wrong reason: any
    # sibling file made .wrapConn() leave the path unwrapped, so the failure
    # was an unrelated SQL parser error about a missing catalog table. Match
    # the message so a future relaxation of that wrap test cannot make an
    # interrupted stream silently readable again.
    expect_error(
        DuckDBDataFrame::DuckDBDataFrame(file.path(tf, "tbl"),
                                         keycol = "__index__"),
        "marked incomplete")
})

test_that("the in-progress marker is written before part 0, not after", {
    # The marker exists to make a crashed write unreadable, so a crash DURING
    # part 0 has to be covered. It previously ran after writeParquet()
    # returned, which left no marker in exactly that case.
    local_mocked_bindings(writeParquet = function(...) stop("simulated crash"),
                          .package = "BiocDuckDB")
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    expect_error(
        writeStreamingResource(
            function(i) if (i == 1L) data.frame(v = 1:3) else NULL,
            path = file.path(tf, "tbl"), dimension = "unbound",
            keycol = NULL, index_max = 3, name = "tbl"),
        "simulated crash")
    expect_true(file.exists(file.path(tf, "tbl", "_INCOMPLETE")))
})

test_that("a completed stream leaves no marker and reads back whole", {
    tf <- tempfile()
    on.exit(unlink(tf, recursive = TRUE), add = TRUE)
    writeStreamingResource(
        function(i) if (i <= 2L) data.frame(v = (1:3) + (i - 1L) * 3L) else NULL,
        path = file.path(tf, "tbl"), dimension = "unbound", keycol = NULL,
        index_max = 6, name = "tbl")

    expect_false(file.exists(file.path(tf, "tbl", "_INCOMPLETE")))
    got <- as.data.frame(
        DuckDBDataFrame::DuckDBDataFrame(file.path(tf, "tbl"),
                                         keycol = "__index__"))
    expect_identical(sort(got$v), 1:6)
})

# --- fixed row-group size (491520) + zonemap stats, via one write --
test_that("coord output uses the fixed 491520 row-group size with stats", {
    n <- 491520L + 10L  # just over one row group
    m <- matrix(seq_len(n), nrow = n, ncol = 1L,
                dimnames = list(paste0("s", seq_len(n)), "f1"))
    d <- tempfile(); dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    writeParquet(m, d)

    pq <- list.files(d, pattern = "\\.parquet$", recursive = TRUE,
                     full.names = TRUE)[1]
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

    groups <- DBI::dbGetQuery(con, sprintf(
        "SELECT DISTINCT row_group_id, row_group_num_rows
           FROM parquet_metadata('%s') ORDER BY row_group_id", pq))
    # multiple bounded groups, the first pinned to the fixed 491520 size
    # (not a fraction of the file) -- so a large resource always parallelizes.
    expect_gte(nrow(groups), 2L)
    expect_identical(as.integer(groups$row_group_num_rows[1]), 491520L)

    # per-row-group min/max statistics are present (the zonemap substrate).
    stats <- DBI::dbGetQuery(con, sprintf(
        "SELECT stats_min, stats_max FROM parquet_metadata('%s')
          WHERE row_group_id = 0 AND stats_min IS NOT NULL", pq))
    expect_gt(nrow(stats), 0L)
})
