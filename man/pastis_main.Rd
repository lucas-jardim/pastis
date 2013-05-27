\name{pastis_main}
\alias{pastis_main}
\title{Phylogenetic Assembly with Soft Taxonomic Inferences}
\usage{
  pastis_main(constraint_tree, taxa_list,
    missing_clades = NA, sequences = NA,
    output_template = NA, output_file = "output.nex",
    paraphyly_constrains = TRUE,
    monophyly_constrains = TRUE, omit_sequences = FALSE)
}
\arguments{
   \item{constraint_tree}{A tree with constraints that are
  forced to be present in all output trees. Either a
  filename to a nexus file readable by read.tree or a ape
  phylo object.}
  
  \item{taxa_list}{A list of all taxa and their clades.
  Either a data frame with columns "taxa" and "clade" or a
  filename for a file readable by read.csv with those
  columns}
  
  \item{missing_clades}{A file containing missing clades.
  Each line of the missing clades file consists of the
  missing clade, the word "include" or "exclude" and a list
  of the reference clades (all separated by commas). Lines
  containing "include" specify that a taxon is contained
  below the MRCA of the reference clades. Lines containing
  "exclude" specify that the missing clade cannot attach
  below the MRCA of the reference clades. #'}
 
  \item{sequences}{A file with all the available sequence
  information in fasta format for details on that format
  see read.dna in the ape package.}

  \item{output_template}{The filename for a template for
  the output nexus file. This file should look like a
  regular mrBayes input file with special tags replacing
  content that will be filled by pastis. In particular:

  <sequences> will be replaced by the sequences (and should
  go below the MATRIX line)

  <ntax> the number of taxa (i.e. "ntax=<ntax>" must be
  somewhere in your template)

  <nchar> the number of characters

  <constraints> the constraints will go here

  <outputfile> where the summaries will be written, (i.e.
  \verb{"sumt filename=<outputfile> burnin ...."} should be
  in your template)

  see default_output_template for an example (which is used
  by default)}
  
   \item{output_file}{The filename for the mrBayes output
  file (will be overwritten if it exists)}

  \item{paraphyly_constrains}{If true, missing clades are
  prevented from entering paraphyletic clades.}

  \item{monophyly_constrains}{If true, missing clades are
  prevented from entering monophyletic clades.}
  
  \item{omit_sequences}{If set to TRUE the sequence file
  (if any) will be ignored. This is useful for testing the
  constraints created by pastis as mrBayes runs much
  quicker without sequence data!}

}
\value{
  NULL
}
\description{
  This function assimilates sequences, taxonomic
  information and tree constraints into a mrBayes file as
  per xxx. This permits the construction of trees that are
  compatible with all of these sources of tdata and contain
  all known taxa.
}
\details{
  This is the main function in pastis which assimilates
  sequences, taxonomic information and tree constraints and
  creates a mrBayes input file. This input file contains
  the tree structure specified by constraint_tree with
  missing taxa in taxa_list and missing clades in
  missing_clades added and placed loosely in the tree using
  the constraint logic outlined in xxx .

  See read_input for a description of the required format
  of the input files. At a minimum the constraining input
  tree and taxa list must be provided.

  In addition to the input checks conducted by read_input
  this function also checks for compatibility between the
  missing genus constraints and constraint tree.

  \preformatted{ }

  \strong{PASTIS}: Phylogenetic Assembly with Soft
  Taxonomic InferenceS? \preformatted{ A bright motmot was
  acting quite rowdy- weaving and squawking quite loudly
  "Pastis is delise" he burped with a sneese "but why is
  everything suddenly cloudy?"} -- \emph{Arne Mooers}
}
\seealso{
  pastis_simple provides a simplified interface to
  pastis_main.

  read_input describes the required file formats

  default_output_template provides an example of the output
  template (also the default)
}

