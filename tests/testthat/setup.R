# Airway counts dataset
data(airway, package = "airway")
airway_counts <- SummarizedExperiment::assay(airway, "counts")
airway_counts_path <- file.path(tempfile(), "airway_counts")
writeParquet(airway_counts, airway_counts_path)


# Helper functions
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
