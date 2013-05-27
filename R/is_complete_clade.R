is_complete_clade <- function (set, clades) 
{
    present_clades <- get_clades(set, clades)
    required_taxa <- unlist(clades[present_clades])
    return(length(intersect(required_taxa, set)) == length(required_taxa))
}