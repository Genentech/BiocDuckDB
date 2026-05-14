### Internal helper functions for inferring Frictionless Data type metadata
### from R/Bioconductor objects. These functions are used by writeParquet()
### to generate complete type information in datapackage.json schemas.
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Field Type Inference
###

.fieldtype <- function(x) {
    UseMethod(".fieldtype")
}

.fieldtype.default <- function(x) {
    list(type = "any")
}

.fieldtype.logical <- function(x) {
    list(type = "boolean")
}

.fieldtype.integer <- function(x) {
    list(type = "integer")
}

.fieldtype.integer64 <- function(x) {
    list(type = "integer", format = "int64")
}

.fieldtype.numeric <- function(x) {
    list(type = "number")
}

.fieldtype.character <- function(x) {
    list(type = "string")
}

.fieldtype.factor <- function(x) {
    list(type = "string",
         categories = I(levels(x)),
         categoriesOrdered = is.ordered(x))
}

.fieldtype.ordered <- function(x) {
    .fieldtype.factor(x)
}

.fieldtype.Date <- function(x) {
    list(type = "date", format = "default") # ISO 8601: YYYY-MM-DD
}

.fieldtype.POSIXct <- function(x) {
    list(type = "datetime", format = "default") # ISO 8601 with timezone
}

.fieldtype.POSIXlt <- function(x) {
    .fieldtype.POSIXct(as.POSIXct(x))
}

.fieldtype.difftime <- function(x) {
    list(type = "duration", format = "default") # ISO 8601 duration
}

.fieldtype.raw <- function(x) {
    list(type = "string", format = "binary") # Base64 encoding in Frictionless
}

.fieldtype.sfc <- function(x) {
    list(type = "string", format = "binary") # WKB geometry
}

.fieldtype.list <- function(x) {
    type <- ifelse(length(x), .fieldtype(x[[1L]])[["type"]], "any")

    field <- list(type = "array",
                  format = "variable",
                  arrayItem = list(type = type))

    lengths_x <- as.vector(unique(lengths(x)))
    if (length(lengths_x) == 1L) {
        field[["format"]] <- "fixed"
        field[["constraints"]] <- list(minLength = lengths_x,
                                       maxLength = lengths_x)
    }

    field
}

.fieldtype.AsIs <- function(x) {
    .fieldtype(unclass(x))
}

.fieldtype.Array <- function(x) {
    arrow_type <- x$type
    type <- arrow_type$ToString()

    type <- if (grepl("^(u?int)", type)) {
        "integer"
    } else if (grepl("^(float|double)", type)) {
        "number"
    } else if (grepl("^(utf8|string|large_)", type)) {
        "string"
    } else if (type == "bool") {
        "boolean"
    } else if (grepl("^date", type)) {
        "date"
    } else if (grepl("^timestamp", type)) {
        "datetime"
    } else if (grepl("^duration", type)) {
        "duration"
    } else if (grepl("^(binary|large_binary)", type)) {
        return(list(type = "string", format = "binary"))
    } else {
        "any"
    }

    result <- list(type = type)

    format <- .arrowTypeToFormat(arrow_type)
    if (!is.null(format)) {
        result[["format"]] <- format
    }

    result
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Arrow Type Selection
###

#' @importFrom arrow infer_type
.arrowType <- function(x) {
    if (is.integer(x)) {
        x <- x[!is.na(x)]
    }

    if (is.integer(x) && length(x) > 0L) {
        DuckDBArray:::.arrowIntType(range(x))
    } else {
        infer_type(x)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Frictionless Format from Arrow Type
###

.arrowTypeToFormat <- function(arrow_type) {
    type_name <- arrow_type$ToString()

    formats <- c("int8" = "int8",
                 "int16" = "int16",
                 "int32" = "int32",
                 "int64" = "int64",
                 "uint8" = "uint8",
                 "uint16" = "uint16",
                 "uint32" = "uint32",
                 "uint64" = "uint64",
                 "binary" = "binary",
                 "large_binary" = "binary")

    if (type_name %in% names(formats)) {
        formats[[type_name]]
    } else {
        NULL
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constraint Helpers
###

.addConstraints <-
function(field, role = NULL, min_value = NULL, max_value = NULL,
         enum_values = NULL, required = NULL, unique = NULL, pattern = NULL)
{
    constraints <- list()

    if (!is.null(role)) {
        constraints <- switch(role,
            "start" = ,
            "end" = list(required = TRUE, minimum = 1L),
            "strand" = list(required = TRUE, enum = c("+", "-", "*")),
            "seqname" = list(required = TRUE),
            "percent" = list(minimum = 0, maximum = 100),
            "count" = list(minimum = 0),
            "index" = list(minimum = 0),
            list()
        )
    }

    if (!is.null(min_value)) constraints[["minimum"]] <- min_value
    if (!is.null(max_value)) constraints[["maximum"]] <- max_value
    if (!is.null(enum_values)) constraints[["enum"]] <- enum_values
    if (!is.null(required)) constraints[["required"]] <- required
    if (!is.null(unique)) constraints[["unique"]] <- unique
    if (!is.null(pattern)) constraints[["pattern"]] <- pattern

    if (length(constraints) > 0L) {
        field[["constraints"]] <- constraints
    }

    field
}

.buildFieldSpec <-
function(name, x, arrow_type = NULL, description = NULL, role = NULL, ...) {
    field <- list(name = name)

    type_spec <- .fieldtype(x)
    field <- c(field, type_spec)

    if (!is.null(arrow_type)) {
        format_str <- .arrowTypeToFormat(arrow_type)
        if (!is.null(format_str) && is.null(field[["format"]])) {
            field[["format"]] <- format_str
        }
    }

    if (!is.null(description)) {
        field[["description"]] <- description
    }

    field <- .addConstraints(field, role = role, ...)

    field
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Schema Enhancement Helpers
###

.addGenomicMetadata <-
function(schema, seqname = "seqnames", start = "start", end = "end",
         strand = "strand")
{
    schema[["genomicCoords"]] <- list(
        seqname = seqname,
        start = start,
        end = end,
        strand = strand
    )

    role_map <- list(
        seqname = "seqname",
        start = "start",
        end = "end",
        strand = "strand"
    )

    for (i in seq_along(schema[["fields"]])) {
        field_name <- schema[["fields"]][[i]][["name"]]
        for (role_name in names(role_map)) {
            if (!is.null(schema[["genomicCoords"]][[role_name]]) &&
                field_name == schema[["genomicCoords"]][[role_name]]) {
                schema[["fields"]][[i]] <- .addConstraints(
                    schema[["fields"]][[i]],
                    role = role_map[[role_name]]
                )
                break
            }
        }
    }

    schema
}

.addGraphMetadata <-
function(schema, nnode, from = "from", to = "to")
{
    # Store nnode as graph metadata
    schema[["graphEdges"]] <- list(
        from = from,
        to = to,
        nnode = nnode
    )

    # Mark from/to columns with semantic roles
    role_map <- list(
        from = "query_node",
        to = "subject_node"
    )

    for (i in seq_along(schema[["fields"]])) {
        field_name <- schema[["fields"]][[i]][["name"]]
        for (role_name in names(role_map)) {
            if (!is.null(schema[["graphEdges"]][[role_name]]) &&
                field_name == schema[["graphEdges"]][[role_name]]) {
                schema[["fields"]][[i]] <- .addConstraints(
                    schema[["fields"]][[i]],
                    role = role_map[[role_name]]
                )
                break
            }
        }
    }

    schema
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Schema Validation Helpers
###

.extractIntType <- function(field) {
    if (field[["type"]] != "integer") {
        return(list(is_int = FALSE))
    }

    format_str <- field[["format"]] %||% "int32"

    unsigned <- grepl("^uint", format_str)
    width <- as.integer(gsub("^u?int", "", format_str))
    constraints <- field[["constraints"]] %||% list()

    list(is_int = TRUE,
         width = width,
         unsigned = unsigned,
         min_value = constraints[["minimum"]],
         max_value = constraints[["maximum"]])
}

.fieldToDuckDBCast <- function(field) {
    field_name <- field[["name"]]
    field_type <- field[["type"]]

    if (field_type == "array" && isTRUE(field[["format"]] == "fixed")) {
        array_size <- field[["constraints"]][["minLength"]]
        if (!is.null(array_size)) {
            item_type <- field[["arrayItem"]][["type"]]
            duckdb_type <- switch(item_type,
                                  "integer" = "INTEGER",
                                  "number" = "DOUBLE",
                                  "string" = "VARCHAR",
                                  "boolean" = "BOOLEAN",
                                  "ANY")
            return(sprintf("%s::%s[%d]", field_name, duckdb_type, array_size))
        }
    }

    if (field_type == "integer") {
        format_str <- field[["format"]]
        if (!is.null(format_str) && format_str != "int32") {
            duckdb_type <- switch(format_str,
                                  "int8" = "TINYINT",
                                  "int16" = "SMALLINT",
                                  "int32" = "INTEGER",
                                  "int64" = "BIGINT",
                                  "uint8" = "UTINYINT",
                                  "uint16" = "USMALLINT",
                                  "uint32" = "UINTEGER",
                                  "uint64" = "UBIGINT",
                                  NULL)
            if (!is.null(duckdb_type)) {
                return(sprintf("%s::%s", field_name, duckdb_type))
            }
        }
    }

    NULL
}
