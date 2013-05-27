create_constraint <- function (type, include, exclude = c(), name = NA, name_class = NA) 
{
    if (length(include) == 1) {
        return()
    }
    type <- tolower(type)
    if (!is.element(type, c("hard", "negative", "partial"))) {
        stop("invalid constraint type specified (valid types are hard, negative and partial)")
    }
    if (!is.na(name_class)) {
        if (!is.element(name_class, names(e$constraint_class_id))) {
            e$constraint_class_id[name_class] <- 0
        }
        e$constraint_class_id[[name_class]] <- e$constraint_class_id[[name_class]] + 
            1
        name <- paste(name_class, e$constraint_class_id[[name_class]], 
            sep = "_")
    }
    if (is.na(name)) {
        name <- as.character(e$constraint_id)
        e$constraint_id <- e$constraint_id + 1
    }
    out <- paste("constraint", name)
    out <- paste(c(out, type, "=", sort(include)), collapse = " ")
    if (length(exclude) > 0) {
        out <- paste(c(out, ":", sort(exclude)), collapse = " ")
    }
    return(out)
}