read_input <- function (constraint_tree, taxa_list, missing_clades = NA, sequences = NA, 
    output_template = NA) 
{
    if (is.na((output_template))) {
        output_template <- default_output_template()
    }
    else {
        if (!file.exists(output_template)) {
            stop(paste("the output template file (", output_template, 
                ") does not exist"))
        }
        output_template <- readChar(output_template, file.info(output_template)$size)
        tags <- c("<ntax>", "<nchar>", "<sequences>", "<constraints>", 
            "<outputfile>")
        for (tag in tags) {
            if (!grep(tag, output_template)) {
                stop(paste("the output template file is missing the tag: ", 
                  tag))
            }
        }
    }
    if (class(constraint_tree) == "character") {
        constraint_tree = read.tree(constraint_tree)
    }
    if (class(constraint_tree) != "phylo") {
        stop("invalid constraint tree specified")
    }
    constraint_tree <- reorder(constraint_tree, "pruningwise")
    constraint_tree <- add_tip_sets(constraint_tree)
    if (class(taxa_list) == "character") {
        taxa_list = read.csv(taxa_list)
    }
    if (class(taxa_list) != "data.frame") {
        stop("invalid taxa list specified")
    }
    if (length(intersect(names(taxa_list), c("taxon", "clade"))) != 
        2) {
        stop("invalid column names in taxa list")
    }
    taxa_list$taxa_char <- as.character(taxa_list$taxon)
    if (length(intersect(constraint_tree$tip.label, taxa_list$taxon)) != 
        length(constraint_tree$tip.label)) {
        missing_taxa = setdiff(constraint_tree$tip.label, intersect(constraint_tree$tip.label, 
            taxa_list$taxon))
        stop(paste("the taxa:", missing_taxa, "is/are present in the constraint tree but not in the taxa list"))
    }
    constraint_taxa <- constraint_tree$tip.label
    constraint_clades <- unique(taxa_list$clade[is.element(taxa_list$taxon, 
        constraint_taxa)])
    if (!is.na(missing_clades)) {
        missing_clades_raw <- scan(missing_clades, what = "character", 
            sep = "\n")
        missing_clades <- list(include = list(), exclude = list())
        for (row in 1:length(missing_clades_raw)) {
            split <- strsplit(missing_clades_raw[row], split = ",")[[1]]
            this_clade <- split[1]
            if (is.element(this_clade, constraint_clades)) {
                stop(paste("the clade", this_clade, "is present in both the constraint tree and the missing clade file"))
            }
            type <- tolower(split[2])
            cogeners <- split[3:length(split)]
            cogeners <- setdiff(cogeners, this_clade)
            if (!is.element(this_clade, taxa_list$clade)) {
                stop(paste("missing clade", this_clade, "is not present in the taxa list"))
            }
            if (length(intersect(cogeners, taxa_list$clade)) != 
                length(cogeners)) {
                discrepancy <- setdiff(cogeners, intersect(cogeners, 
                  taxa_list$clade))
                stop(paste("missing clade", this_clade, "has sister clade(s) that are not in the taxa list"))
            }
            if (!is.element(type, c("include", "exclude"))) {
                stop(paste("missing clade", this_clade, "has an invalid constraint type (valid types are include and exclude)"))
            }
            missing_cogeners <- setdiff(cogeners, constraint_clades)
            if (length(missing_cogeners) > 0) {
                warning(paste("the missing clade constraint for", 
                  this_clade, "refers to missing clade(s):", 
                  paste(missing_cogeners, collapse = ","), ", this reference will be ignored when positioning this clade"))
            }
            cogeners <- intersect(cogeners, constraint_clades)
            missing_clades[[type]][[this_clade]] <- cogeners
        }
    }
    else {
        missing_clades <- list(include = list(), exclude = list())
    }
    if (is.na(sequences)) {
        sequences <- matrix(c(1, 1))
    }
    else {
        e <- try(sequences <- read.dna(sequences, format = "fasta", 
            as.character = TRUE))
        if (class(e) == "try-error") {
            stop(paste("could not load sequence data from", sequences, 
                ". See read.dna for a description of the fasta format."))
        }
    }
    return(list(constraint_tree = constraint_tree, taxa_list = taxa_list, 
        missing_clades = missing_clades, sequences = sequences, 
        output_template = output_template))
}