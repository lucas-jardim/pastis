get_clades <- function (set, clades) 
{
    out = c()
    for (clade in names(clades)) {
        if (length(intersect(set, clades[[clade]])) > 0) {
            out <- c(out, clade)
        }
    }
    return(out)
}