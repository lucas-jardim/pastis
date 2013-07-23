\name{pastis_simple}
\alias{pastis_simple}
\title{A simplified interface to the main pastis function.}
\usage{
  pastis_simple(pastisData = NULL, base_name,
    paraphyly_constrains = TRUE,
    monophyly_constrains = TRUE, omit_sequences = FALSE)
}
\arguments{
  \item{base_name}{The base name for all input files may
  include a leading directory, but should not include a
  trailing .}

  \item{pastisData}{Input data object of class pastisData}

  \item{paraphyly_constrains}{If TRUE, missing clades are
  prevented from entering paraphyletic clades.}

  \item{monophyly_constrains}{If TRUE, missing clades are
  precented from entering monophyletic clades.}

  \item{omit_sequences}{If set to TRUE the sequence file
  (if any) will be ignored. This is useful for testing the
  constraints created by pastis as mrBayes runs much
  quicker without sequence data!}
}
\description{
  This function assimilates sequences, taxonomic
  information and tree constraints into a mrBayes file as
  per xxx. This permits the construction of trees that are
  compatible with all of these sources of tdata and contain
  all known taxa.
}
\details{
  This is a simplified version of pastis_main which assumes
  that (i) the data are input as a pastisData object (e.g.
  created with read_input), or (ii) the input files all
  have the same base (specified by base_name) with
  diffferent extensions for each file type:

  .sequences: the sequence file (FASTA format)

  .tree: the constraint tree (Newick format)

  .taxa: a list of taxa (see below)

  .missingclades: locations of missing clades (see below)

  .template: the template for the mrBayes input file (see
  default_output_template() )

  The taxa file consists of a header 'taxon,clade' with
  each subsequent line containing a taxon,clade pair
  (separated by a comma)

  Each line of the missingclades file consists of the
  missing clade, the word include or exclude and a list of
  the reference clades (all separated by commas). Lines
  containing include specify that a taxon is contained
  below the MRCA of the reference clades. Lines containing
  exclude specify that the missing clade cannot attach
  below the MRCA of the reference clades.

  The mrBayes input file is written to base_name.nexus.

  \preformatted{ }

  \strong{PASTIS}: Phylogenetic Assembly with Soft
  Taxonomic InferenceS? \preformatted{ A bright motmot was
  acting quite rowdy- weaving and squawking quite loudly
  "Pastis is delise" he burped with a sneese "but why is
  everything suddenly cloudy?"} -- \emph{Arne Mooers}
}
\examples{
# Generate MrBayes input files with constraints
\dontrun{
data(accipitridaeFullPastis)
pastis_simple(accipitridaeFullPastis, base_name="Accipitridae")

data(accipitridaeBasicPastis)
pastis_simple(accipitridaeBasicPastis, base_name="AccipitridaeBasic")
}

data(pastis_data_1)
pastis_simple(pastis_data_1, base_name="pastis_data_1")

data(pastis_data_2)
pastis_simple(pastis_data_2, base_name="pastis_data_2")

data(pastis_data_3)
pastis_simple(pastis_data_3, base_name="pastis_data_3")
}

