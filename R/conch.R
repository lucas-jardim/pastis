conch <- function (constraint_tree, mrbayes_output, simple_edge_scaling = TRUE, 
    species_set = NA) 
{
    if (class(constraint_tree) == "character") {
        constraint_tree = read.tree(constraint_tree)
    }
    if (class(constraint_tree) != "phylo") {
        stop("invalid constraint tree specified")
    }
    constraining_taxa <- constraint_tree$tip.label
    mrbayes_trees <- read.nexus(mrbayes_output)
    missing_taxa <- setdiff(mrbayes_trees[[1]]$tip.label, constraining_taxa)
    if (is.na(species_set)) {
        species_set <- missing_taxa
    }
    count <- 1
    for (taxon in species_set) {
        print(paste("Processing species ", count, " of ", length(species_set), 
            ": ", taxon, sep = ""))
        count <- count + 1
        constraint_tree$edge.length <- c(1:dim(constraint_tree$edge)[1])
        constraint_tree$edge.length[] <- 0
        constraint_tree$root.edge <- 0
        for (tree_index in 1:length(mrbayes_trees)) {
            this_tree <- mrbayes_trees[[tree_index]]
            this_tree <- drop.tip(this_tree, setdiff(this_tree$tip.label, 
                c(constraining_taxa, taxon)))
            edge <- this_tree$edge[, 2] == which(this_tree$tip.label == 
                taxon)
            this_tree <- extract.clade(this_tree, this_tree$edge[this_tree$edge[, 
                2] == which(this_tree$tip.label == taxon), 1])
            mrca_descendants <- intersect(this_tree$tip.label, 
                constraining_taxa)
            if (length(mrca_descendants) == 1) {
                mrca_edge <- which.edge(constraint_tree, mrca_descendants)[1]
            }
            else {
                mrca_node <- constraint_tree$edge[which.edge(constraint_tree, 
                  mrca_descendants)[1], 1]
                mrca_edge <- which(constraint_tree$edge[, 2] == 
                  mrca_node)
            }
            if (simple_edge_scaling) {
                constraint_tree$edge.length[mrca_edge] <- 1
            }
            else {
                constraint_tree$edge.length[mrca_edge] <- constraint_tree$edge.length[mrca_edge] + 
                  1
            }
        }
        write.tree(constraint_tree, file = paste("taxonposition_", 
            taxon, ".tree", sep = ""))
    }
}