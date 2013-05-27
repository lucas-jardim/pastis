find_cherry <- function (tree) 
{
    n_species <- length(tree$tip.label)
    for (edge in n_species + 1:max(tree$edge)) {
        select <- tree$edge[, 1] == edge
        if (sum(tree$edge[select, 2] > n_species) == 0) {
            break
        }
    }
    return(edge)
}