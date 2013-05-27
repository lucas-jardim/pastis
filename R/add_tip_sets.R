add_tip_sets <- function (tree) 
{
    tip.sets <- list()
    for (species in 1:length(tree$tip.label)) {
        tip.sets[[species]] <- tree$tip.label[species]
    }
    tree$tip.sets <- tip.sets
    return(tree)
}