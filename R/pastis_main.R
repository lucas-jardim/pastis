pastis_main <- function (constraint_tree, taxa_list, missing_clades = NA, sequences = NA, 
    output_template = NA, output_file = "output.nex", paraphyly_constrains = TRUE, 
    monophyly_constrains = TRUE, omit_sequences = FALSE) 
{
    constraints <- c()
    reset_constraint_counters()
    input <- read_input(constraint_tree, taxa_list, missing_clades, 
        sequences, output_template)
    constraint_tree <- input$constraint_tree
    taxa_list <- input$taxa_list
    missing_clades <- input$missing_clades
    sequences <- input$sequences
    output_template <- input$output_template
    if (omit_sequences) {
        sequences = NA
    }
    clade_list = list()
    clade_list_constraining = list()
    for (clade in levels(taxa_list$clade)) {
        clade_list[[clade]] <- as.character(taxa_list$taxon[taxa_list$clade == 
            clade])
        clade_list_constraining[[clade]] <- intersect(taxa_list$taxon[taxa_list$clade == 
            clade], constraint_tree$tip.label)
        if (length(clade_list_constraining[[clade]]) == 0 & length(clade_list[[clade]]) > 
            1) {
            constraints <- c(constraints, create_constraint("hard", 
                clade_list[[clade]], name_class = "missing_clade"))
        }
    }
    remaining_clades <- names(clade_list)
    constraining_taxa <- constraint_tree$tip.label
    all_taxa <- as.character(taxa_list$taxon)
    tree <- constraint_tree
    st <- FALSE
    pendant_edge_index <- 0
    while (length(tree$tip.label) > 1) {
        if (pendant_edge_index < length(tree$tip.label)) {
            pendant_edge_index <- pendant_edge_index + 1
            index <- pendant_edge_index
        }
        else {
            edge <- find_cherry(tree)
            n_leaves <- sum(tree$edge[, 1] == edge)
            found <- FALSE
            if (n_leaves > 2) {
                for (s in 2:(n_leaves - 1)) {
                  for (combination in combn(seq(n_leaves), s, 
                    simplify = FALSE)) {
                    collapsed_temp <- collapse_cherry(tree, edge, 
                      combination)
                    cherry_taxa <- collapsed_temp[[1]]$tip.sets[[collapsed_temp[[2]]]]
                    if (is_complete_clade(cherry_taxa, clade_list_constraining)) {
                      collapsed <- collapse_cherry(tree, edge, 
                        combination)
                      tree <- collapsed[[1]]
                      index <- collapsed[[2]]
                      found <- TRUE
                      break
                    }
                  }
                  if (found == TRUE) {
                    break
                  }
                }
            }
            if (found == FALSE) {
                collapsed <- collapse_cherry(tree, edge)
                tree <- collapsed[[1]]
                index <- collapsed[[2]]
            }
        }
        cherry_taxa <- tree$tip.sets[[index]]
        complete <- is_complete_clade(cherry_taxa, clade_list_constraining)
        if (complete) {
            clades <- get_clades(cherry_taxa, clade_list_constraining)
            new_clades <- intersect(clades, remaining_clades)
            remaining_clades <- setdiff(remaining_clades, new_clades)
            cherry_all_taxa <- unique(c(unlist(clade_list[clades]), 
                cherry_taxa))
            includes_missing_taxa <- (length(cherry_taxa) != 
                length(cherry_all_taxa))
            if ((length(clades) == 1 & monophyly_constrains) | 
                (length(new_clades) > 0 & length(clades) > 1 & 
                  paraphyly_constrains)) {
                constraints <- c(constraints, create_constraint("hard", 
                  cherry_all_taxa, name_class = "ctree"))
            }
            else {
                exclude <- setdiff(constraining_taxa, cherry_taxa)
                if (length(exclude) > 0) {
                  constraints <- c(constraints, create_constraint("partial", 
                    cherry_all_taxa, exclude = exclude, name_class = "ctree"))
                }
            }
            for (type in c("include", "exclude")) {
                for (mc in names(missing_clades[[type]])) {
                  if (prod(is.element(missing_clades[[type]][[mc]], 
                    clades)) == 1) {
                    if (type == "exclude") {
                      constraints <- c(constraints, create_constraint("partial", 
                        cherry_all_taxa, exclude = clade_list[[mc]], 
                        name_class = paste("missing_clade_exclude_", 
                          mc, sep = "")))
                      if (!is.element(mc, names(missing_clades[["include"]]))) {
                        stop(paste("Missing clade", mc, "reached exclude statement without a remaining include"))
                      }
                    }
                    else {
                      cherry_all_taxa <- c(cherry_all_taxa, clade_list[[mc]])
                      constraints <- c(constraints, create_constraint("partial", 
                        cherry_all_taxa, exclude = setdiff(constraining_taxa, 
                          cherry_all_taxa), name_class = paste("missing_clade_include_", 
                          mc, sep = "")))
                      if (is.element(mc, names(missing_clades[["exclude"]]))) {
                        stop(paste("Missing clade", mc, "reached include statement while exclude statements remain"))
                      }
                    }
                    missing_clades[[type]][[mc]] <- NULL
                  }
                }
            }
            names(cherry_all_taxa) <- NULL
            tree$tip.sets[[index]] <- cherry_all_taxa
        }
        else {
            constraints <- c(constraints, create_constraint("partial", 
                cherry_taxa, exclude = setdiff(constraining_taxa, 
                  cherry_taxa), name_class = "ctree"))
        }
    }
    constraint_string <- paste(constraints, collapse = ";\n")
    constraint_splits <- strsplit(constraints, split = " ")
    constraint_names <- paste(sapply(constraint_splits, function(x) {
        return(x[2])
    }), collapse = ", ")
    constraint_string <- paste(constraint_string, paste("prset topologypr=constraints(", 
        constraint_names, ")", sep = ""), sep = ";\n")
    constraint_string <- paste(constraint_string, ";", sep = "")
    sequence_string = ""
    n_char <- dim(sequences)[2]
    missing_sequence <- paste(rep("?", n_char), collapse = "")
    for (taxon in all_taxa) {
        if (is.element(taxon, row.names(sequences))) {
            sequence_string <- paste(sequence_string, "\n", taxon, 
                paste(sequences[taxon, ], collapse = ""), sep = " ")
        }
        else {
            sequence_string <- paste(sequence_string, "\n", taxon, 
                " ", missing_sequence, sep = "")
        }
    }
    if (grep("<ntax>", output_template)) {
        output <- output_template
    }
    else {
        output <- readChar(output_template, file.info(output_template)$size)
    }
    output <- sub("<ntax>", length(all_taxa), output)
    output <- sub("<nchar>", n_char, output)
    output <- sub("<outputfile>", output_file, output)
    output <- sub("<constraints>", constraint_string, output)
    output <- sub("<sequences>", sequence_string, output)
    writeChar(output, output_file)
}