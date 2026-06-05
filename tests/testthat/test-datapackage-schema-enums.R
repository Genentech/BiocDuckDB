# DRY guard: the controlled vocabularies in the bundled profile
# (inst/schema/biocduckdb-profile.json) must stay in sync with what writeParquet()
# can emit. The reference sets below are the authoritative vocabularies; if a new
# layout/dimension/model is added to writeParquet() (or removed), update both the
# profile enum and the matching reference set here, and this test keeps them
# locked together. Cross-checked against the layout/dimension defaults and
# package=list(model=) values in R/writeParquet.R.
# Run: library(BiocDuckDB); library(testthat); source("test-datapackage-schema-enums.R")

# Resource layouts writeParquet() can record (one per primitive / experiment method).
.KNOWN_LAYOUTS <- c(
    "coord_array",
    "data_frame",
    "transposed_data_frame",
    "embedding_table",
    "genomic_ranges",
    "genomic_ranges_list",
    "graph_edges",
    "spatial_points",
    "spatial_shapes",
    "nested_data_frame",
    "nested_experiment"
)

# Biological axes (match.arg dimension defaults across the methods).
.KNOWN_DIMENSIONS <- c("feature", "sample", "crossed", "unbound")

# Package-level container models (package=list(model=) plus the RangedSE ifelse;
# experiment_list is built internally and writes no top-level datapackage.json).
.KNOWN_MODELS <- c(
    "summarized_experiment",
    "ranged_summarized_experiment",
    "single_cell_experiment",
    "experiment_list",
    "multi_assay_experiment",
    "multi_assay_spatial_experiment"
)

.profileEnums <- function() {
    skip_if_not_installed("jsonlite")
    profile <- system.file("schema", "biocduckdb-profile.json", package = "BiocDuckDB")
    if (!nzchar(profile) || !file.exists(profile)) {
        skip("bundled schema profile not found (package not installed with inst/schema)")
    }
    spec <- jsonlite::fromJSON(profile, simplifyVector = FALSE)
    resource <- spec[["$defs"]][["biocduckdb_resource"]][["properties"]]
    list(
        model = unlist(spec[["properties"]][["model"]][["enum"]]),
        layout = unlist(resource[["layout"]][["enum"]]),
        dimension = unlist(resource[["dimension"]][["enum"]])
    )
}

test_that("profile layout enum matches the writeParquet() layout vocabulary", {
    enums <- .profileEnums()
    expect_setequal(enums$layout, .KNOWN_LAYOUTS)
})

test_that("profile dimension enum matches the biological axes", {
    enums <- .profileEnums()
    expect_setequal(enums$dimension, .KNOWN_DIMENSIONS)
})

test_that("profile model enum matches the documented container models", {
    enums <- .profileEnums()
    expect_setequal(enums$model, .KNOWN_MODELS)
})
