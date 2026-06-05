# DRY guard: the CLOSED vocabularies in the bundled profile
# (inst/schema/biocduckdb-profile.json) must stay in sync with what writeParquet()
# can emit. 'layout' and 'dimension' are closed enums and are checked setequal
# against the reference sets below; if a new value is added to writeParquet(),
# update both the profile enum and the matching reference set here. 'model' is an
# OPEN string (general galaxy-schema framework), so it is checked for openness and
# its documented reader values, not pinned to an enum.
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

# Reader values BiocDuckDB ships dedicated container methods for. 'model' is an
# OPEN string (BiocDuckDB is a general galaxy-schema framework: any reconstruction
# key is permitted, unknown/absent yields a SimpleList), so these are documented
# in the profile description, NOT enforced as an enum.
.DOCUMENTED_MODELS <- c(
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
        model = spec[["properties"]][["model"]],
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

test_that("profile model is an open string documenting the known readers", {
    enums <- .profileEnums()
    # model must NOT be a closed enum (general galaxy-schema framework).
    expect_identical(enums$model[["type"]], "string")
    expect_null(enums$model[["enum"]])
    # The known reader values should be documented in the description so the
    # vocabulary stays discoverable even though it is not enforced.
    desc <- enums$model[["description"]]
    expect_true(is.character(desc) && nzchar(desc))
    for (m in .DOCUMENTED_MODELS) {
        expect_true(grepl(m, desc, fixed = TRUE),
                    info = paste("model description should mention", m))
    }
})
