# Airway counts dataset
data(airway, package = "airway")
airway_counts <- SummarizedExperiment::assay(airway, "counts")
airway_counts_path <- file.path(tempfile(), "airway_counts")
BiocDuckDB::writeParquet(airway_counts, airway_counts_path)


# Helper functions
checkDuckDBTable <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBTable")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_gte(nrow(object), nrow(expected))
    expect_gte(NROW(object), NROW(expected))
    expect_equal(nkey(object) + ncol(object), ncol(expected))
    expect_equal(nkey(object) + NCOL(object), NCOL(expected))
    expected_cols <- c(setdiff(colnames(expected), keynames(object)), keynames(object))
    expect_identical(c(colnames(object), keynames(object)), expected_cols)
    if (nkey(object) == 0L) {
        object <- as.data.frame(object)
        expect_gte(nrow(object), nrow(expected))
        expect_equal(ncol(object) - 1L, ncol(expected))
    } else {
        df <- as.data.frame(object)
        df <- df[match(do.call(paste, expected[, keynames(object), drop = FALSE]),
                       do.call(paste, df[, keynames(object), drop = FALSE])), ]
        rownames(df) <- NULL
        dcol_names <- setdiff(colnames(expected), keynames(object))
        expected <- expected[, c(dcol_names, keynames(object)), drop = FALSE]
        expect_equivalent(df, expected)
    }
}

checkDuckDBDataFrame <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBDataFrame")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(nrow(object), nrow(expected))
    expect_identical(ncol(object), ncol(expected))
    expect_setequal(rownames(object), rownames(expected))
    expect_identical(colnames(object), colnames(expected))
    if (nkey(object) == 0L) {
        object <- as.data.frame(object)
        expect_identical(nrow(object), nrow(expected))
        expect_identical(ncol(object), ncol(expected))
        expect_identical(colnames(object), colnames(expected))
    } else {
        expect_identical(as.data.frame(object)[rownames(expected), , drop=FALSE], expected)
    }
}

checkDuckDBColumn <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBColumn")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(length(object), length(expected))
    if (nkey(object@table) == 0L) {
        object <- as.vector(object)
        expect_identical(length(object), length(expected))
    } else {
        expect_identical(names(object), names(expected))
        expect_equal(as.vector(object), expected)
        expect_equal(realize(object), expected)
    }
}

checkDuckDBAtomicList <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBAtomicList")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(length(object), length(expected))
    expect_identical(elementNROWS(object), elementNROWS(expected))
    if (nkey(object@table) == 0L) {
        object <- as.list(object)
        expect_identical(length(object), length(expected))
    } else {
        expect_identical(names(object), names(expected))
        expect_equal(as.list(object), expected)
        expect_equal(realize(object), expected)
    }
}

checkDuckDBEmbeddings <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBEmbeddings")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(dim(object), dim(expected))
    expect_identical(nrow(object), nrow(expected))
    expect_identical(ncol(object), ncol(expected))
    if (nkey(object@table) == 0L) {
        mat <- as.matrix(object)
        expect_identical(dim(mat), dim(expected))
    } else {
        expect_identical(rownames(object), rownames(expected))
        expect_equal(as.matrix(object), expected, tolerance = 1e-10)
    }
}

checkDuckDBTransposedDataFrame <- function(object, texpected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBTransposedDataFrame")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(nrow(object), ncol(texpected))
    expect_identical(ncol(object), nrow(texpected))
    expect_identical(rownames(object), colnames(texpected))
    expect_setequal(colnames(object), rownames(texpected))
    if (nkey(t(object)) == 0L) {
        tobject <- as.data.frame(t(object))
        expect_identical(nrow(tobject), nrow(texpected))
        expect_identical(ncol(tobject), ncol(texpected))
        expect_identical(colnames(tobject), colnames(texpected))
    } else {
        expect_identical(as.data.frame(t(object))[rownames(texpected), , drop=FALSE], texpected)
    }
}

checkDuckDBDataFrameList <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBDataFrameList")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(length(object), length(expected))
    expect_identical(names(object), names(expected))
    expect_identical(NROW(object), NROW(expected))
    expect_identical(ROWNAMES(object), ROWNAMES(expected))
    expect_identical(elementNROWS(object), elementNROWS(expected))
    expect_identical(nrows(object), nrows(expected))
    expect_identical(ncols(object), ncols(expected))
    expect_identical(dims(object), dims(expected))
    for (i in seq_along(object)) {
        expect_setequal(rownames(object)[[i]], rownames(expected)[[i]])
    }
    expect_identical(colnames(object), colnames(expected))
    expect_identical(mcols(object), mcols(expected))
    expect_identical(columnMetadata(object), columnMetadata(expected))
    expect_identical(commonColnames(object), commonColnames(expected))
    checkDuckDBDataFrame(unlist(object), as.data.frame(unlist(expected, use.names = FALSE)))
}

checkDuckDBSelfHits <- function(object, expected) {
    expect_s4_class(object, "DuckDBSelfHits")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(length(object), length(expected))
    expect_identical(nnode(object), nnode(expected))
    expect_identical(nLnode(object), nLnode(expected))
    expect_identical(nRnode(object), nRnode(expected))
    expect_identical(countLnodeHits(object), countLnodeHits(expected))
    expect_identical(countRnodeHits(object), countRnodeHits(expected))
    if (length(object) > 0L) {
        expect_identical(dbconn(object), acquireDuckDBConn())
        expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    }
    if (nkey(object@frame) > 0L) {
        expect_identical(unname(as.integer(from(object))), as.integer(from(expected)))
        expect_identical(unname(as.integer(to(object))), as.integer(to(expected)))
        expect_identical(as(object, "SelfHits"), expected)
        DF <- as(object, "DFrame")
        rownames(DF) <- NULL
        expect_identical(DF, as(expected, "DFrame"))
        df <- as.data.frame(expected)
        rownames(df) <- NULL
        expect_identical(df, as.data.frame(expected))
    }
}

checkDuckDBMatrix <- function(object, expected) {
    expect_true(validObject(object))
    expect_s4_class(object, "DuckDBMatrix")
    expect_identical(dbconn(object), acquireDuckDBConn())
    expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    expect_identical(type(object), typeof(expected))
    expect_identical(length(object), length(expected))
    expect_identical(dim(object), dim(expected))
    expect_identical(dimnames(object), dimnames(expected))
    expect_equal(as.matrix(object), expected)
    expect_equal(as(object, "CsparseMatrix"), as(expected, "CsparseMatrix"))
    expect_equivalent(as(object, "SparseMatrix"), as(expected, "COO_SparseMatrix"))
    expect_equivalent(as(object, "COO_SparseMatrix"), as(expected, "COO_SparseMatrix"))
}

checkDuckDBGRanges <- function(object, expected) {
    expect_s4_class(object, "DuckDBGRanges")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(length(object), length(expected))
    if (length(object) > 0L) {
        expect_identical(dbconn(object), acquireDuckDBConn())
        expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    }
    if (nkey(object@frame) > 0L) {
        expect_setequal(names(object), names(expected))
        object <- object[names(expected)]
        expect_identical(unname(as.vector(seqnames(object))), as.character(seqnames(expected)))
        expect_identical(unname(as.vector(start(object))), start(expected))
        expect_identical(unname(as.vector(end(object))), end(expected))
        expect_identical(unname(as.vector(width(object))), width(expected))
        expect_identical(unname(as.vector(strand(object))), as.character(strand(expected)))
        expect_setequal(seqlevels(object), seqlevels(expected))
        expect_setequal(seqlengths(object), seqlengths(expected))
        expect_setequal(isCircular(object), isCircular(expected))
        expect_setequal(genome(object), genome(expected))
        df <- as.data.frame(expected)
        for (j in names(df)) {
            if (is.factor(df[[j]])) {
                df[[j]] <- as.character(df[[j]])
            }
        }
        expect_identical(as.data.frame(object)[names(expected), , drop=FALSE], df)
    }
}

checkDuckDBGRangesList <- function(object, expected) {
    expect_s4_class(object, "DuckDBGRangesList")
    expect_true(length(capture.output(show(object))) > 0L)
    expect_identical(length(object), length(expected))
    expect_identical(names(object), names(expected))
    expect_identical(elementNROWS(object), elementNROWS(expected))
    if (length(object) > 0L) {
        expect_identical(dbconn(object), acquireDuckDBConn())
        expect_s3_class(tblconn(object), "tbl_duckdb_connection")
    }
    if (nkey(object@frame) > 0L) {
        expect_identical(as.list(seqnames(object)), lapply(as.list(seqnames(expected)), as.character))
        expect_identical(as.list(start(object)), as.list(start(expected)))
        expect_identical(as.list(end(object)), as.list(end(expected)))
        expect_identical(as.list(width(object)), as.list(width(expected)))
        expect_identical(as.list(strand(object)), lapply(as.list(strand(expected)), as.character))
        expect_setequal(seqlevels(object), seqlevels(expected))
        expect_setequal(seqlengths(object), seqlengths(expected))
        expect_setequal(isCircular(object), isCircular(expected))
        expect_identical(as(mcols(object), "DFrame"), mcols(expected))
    }
}

checkDuckDBDataFrameSpatial <- function(object, expected, coords = c("x", "y")) {
    checkDuckDBDataFrame(object, expected)
    for (cn in coords) {
        if (cn %in% colnames(expected))
            expect_true(cn %in% colnames(object))
    }
}

makeSpatialMASEFixture <- function() {
    skip_if_not_installed("MultiAssaySpatialExperiment")
    skip_if_not_installed("sf")
    pts <- S4Vectors::DataFrame(
        x = c(1.5, 2.5, 3.5),
        y = c(1.5, 2.5, 3.5),
        instance_id = c("A", "B", "C"))
    shp_df <- S4Vectors::DataFrame(
        instance_id = c("cell1", "cell2", "cell3"),
        geometry = sf::st_sfc(
            sf::st_polygon(list(matrix(c(1,1,2,1,2,2,1,2,1,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,1,3,1,3,2,2,2,2,1), ncol=2, byrow=TRUE))),
            sf::st_polygon(list(matrix(c(2,2,3,2,3,3,2,3,2,2), ncol=2, byrow=TRUE)))))
    MultiAssaySpatialExperiment::MultiAssaySpatialExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            assay1 = matrix(c(1:9), 3, 3,
                dimnames = list(paste0("G", 1:3), c("A", "B", "C")))),
        colData = S4Vectors::DataFrame(row.names = "s1"),
        sampleMap = S4Vectors::DataFrame(
            assay = "assay1", primary = "s1", colname = c("A", "B", "C")),
        points = MultiAssaySpatialExperiment::PointsLayerList(centroids = pts),
        shapes = MultiAssaySpatialExperiment::ShapesLayerList(cells = shp_df),
        spatialMap = S4Vectors::DataFrame(
            assay = "assay1", colname = c("A", "B", "C"),
            element_type = "points", region = "centroids",
            instance_id = c("A", "B", "C")),
        metadata = list(transforms = list(
            "points/centroids" = list(
                global = list(type = "identity"),
                scaled = list(type = "scale", scale = c(2, 2))),
            "shapes/cells" = list(
                global = list(type = "identity")))))
}

makeLazySpatialMASE <- function(path = NULL) {
    skip_if_not_installed("MultiAssaySpatialExperiment")
    mase <- makeSpatialMASEFixture()
    if (is.null(path))
        path <- file.path(tempdir(), paste0("mase_lazy_", sample.int(1e6, 1L)))
    BiocDuckDB::writeParquet(mase, path)
    BiocDuckDB::readParquet(path)
}
