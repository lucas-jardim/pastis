collapse_cherry <- function (tree, edge = NA, combine = NA) 
{
    if (is.na(edge)) {
        n_species <- length(tree$tip.label)
        edge <- find_cherry(tree)
        if (!is.na(combine)) {
            stop("You can only specify which branches to combine if you have specified an edge")
        }
    }
    if (!is.element("tip.sets", names(tree))) {
        tree <- add_tip_sets(tree)
    }
    rows <- tree$edge[, 1] == edge
    if (!is.na(combine)[1]) {
        if (max(combine) > sum(rows) | min(combine) < 1) {
            stop("Invalid edges specified in combine")
        }
        if (length(unique(combine)) != length(combine)) {
            stop("Elements in combine must be unique")
        }
        if (length(combine) < 2) {
            stop("You must specify multiple edges to combine")
        }
    }
    if (length(combine) == sum(rows) | is.na(combine)[1]) {
        leaves <- sort(tree$edge[rows, 2])
        new_leaf <- leaves[1]
        new_edge <- tree$edge
        new_edge[new_edge[, 2] == edge, 2] = new_leaf
        new_edge <- new_edge[!rows, ]
        select <- new_edge[, 2] >= edge
        new_edge[select, 2] <- new_edge[select, 2] - 1
        select <- new_edge[, 1] >= edge
        new_edge[select, 1] <- new_edge[select, 1] - 1
        for (removed_leaf in seq(length(leaves), 2, by = -1)) {
            select <- new_edge[, 2] >= leaves[removed_leaf]
            new_edge[select, 2] <- new_edge[select, 2] - 1
            select <- new_edge[, 1] >= leaves[removed_leaf]
            new_edge[select, 1] <- new_edge[select, 1] - 1
        }
        new_tip.label <- tree$tip.label[-leaves[2:length(leaves)]]
        new_tip.label[leaves[1]] <- "merge"
        new_tip.sets <- tree$tip.sets
        new_tip.sets[[leaves[1]]] <- unique(unlist(new_tip.sets[leaves]))
        new_tip.sets <- new_tip.sets[-leaves[2:length(leaves)]]
        new_tree <- tree
        new_tree$tip.label <- new_tip.label
        new_tree$edge <- new_edge
        new_tree$Nnode <- new_tree$Nnode - 1
        new_tree$tip.sets <- new_tip.sets
    }
    else {
        rows <- which(tree$edge[, 1] == edge)[combine]
        leaves <- sort(tree$edge[rows, 2])
        new_edge <- tree$edge
        new_edge <- new_edge[-rows[-1], ]
        for (removed_leaf in seq(length(leaves), 2, by = -1)) {
            select <- new_edge[, 2] >= leaves[removed_leaf]
            new_edge[select, 2] <- new_edge[select, 2] - 1
            select <- new_edge[, 1] >= leaves[removed_leaf]
            new_edge[select, 1] <- new_edge[select, 1] - 1
        }
        new_tip.label <- tree$tip.label[-leaves[2:length(leaves)]]
        new_tip.label[leaves[1]] <- "merge"
        new_tip.sets <- tree$tip.sets
        new_tip.sets[[leaves[1]]] <- unique(unlist(new_tip.sets[leaves]))
        new_tip.sets <- new_tip.sets[-leaves[2:length(leaves)]]
        new_tree <- tree
        new_tree$tip.label <- new_tip.label
        new_tree$edge <- new_edge
        new_tree$Nnode <- new_tree$Nnode
        new_tree$tip.sets <- new_tip.sets
    }
    return(list(new_tree, leaves[1]))
}