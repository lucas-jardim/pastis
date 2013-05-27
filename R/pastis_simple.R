pastis_simple <- function (base_name, paraphyly_constrains = TRUE, monophyly_constrains = TRUE, 
    omit_sequences = FALSE) 
{
    file_check <- function(base, extension, silent = FALSE) {
        file <- paste(base, extension, sep = ".")
        if (file.exists(file)) {
            return(file)
        }
        else {
            if (!silent) {
                warning(paste("No", extension, "file was found, proceeding without (this is fine if it is intentional)"))
            }
            return(NA)
        }
    }
    return(pastis_main(constraint_tree = file_check(base_name, 
        "tree"), taxa_list = file_check(base_name, "taxa"), missing_clades = file_check(base_name, 
        "missingclades"), sequences = file_check(base_name, "sequences"), 
        output_file = paste(base_name, "nexus", sep = "."), output_template = file_check(base_name, 
            "template"), paraphyly_constrains = paraphyly_constrains, 
        monophyly_constrains = monophyly_constrains, omit_sequences = omit_sequences))
}