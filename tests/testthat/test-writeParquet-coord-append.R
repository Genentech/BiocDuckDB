test_that("writeParquet forwards append arguments to writeCoordArray", {
    path <- tempfile()
    dir.create(path)

    m1 <- matrix(1:20, nrow = 4L, ncol = 5L)
    grid1 <- S4Arrays::RegularArrayGrid(dim(m1), c(2L, 2L))
    writeParquet(m1, path, grid = grid1, grid_suffix = "group__")

    pq_before <- list.files(path, pattern = "\\.parquet$", recursive = TRUE)

    m2 <- matrix(100:119, nrow = 4L, ncol = 5L)
    grid2 <- S4Arrays::RegularArrayGrid(dim(m2), c(2L, 2L))
    suppressWarnings(writeParquet(m2, path,
                                  grid = grid2,
                                  grid_suffix = "group__",
                                  append = TRUE,
                                  along = 2L,
                                  offset = ncol(m1),
                                  group_offset = dim(grid1)[2L]))

    pq_after <- list.files(path, pattern = "\\.parquet$", recursive = TRUE)
    expect_gt(length(pq_after), length(pq_before))
    expect_true(any(grepl("index2group__=3", pq_after)))
})
