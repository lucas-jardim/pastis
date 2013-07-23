\name{conch}
\alias{conch}
\title{mrBayes output interrogator (COnstraint CHecker)}
\usage{
  conch(constraint_tree, mrbayes_output,
    simple_edge_scaling = TRUE, species_set = NA)
}
\arguments{
  \item{constraint_tree}{the constraint tree used with
  pastis (either the filename or ape phylo tree)}

  \item{mrbayes_output}{the mrBayes output .t file (either
  the filename or ape multiPhylo tree object)}

  \item{simple_edge_scaling}{boolean, see function
  description}

  \item{species_set}{if specified, this should be a list of
  species and output trees will only be generated for these
  species (the default is all missing species)}
}
\value{
  NULL. A file for each missing taxon is created in the
  current directory
}
\description{
  This function examines mrBayes output from a pastis run
  to examine where in the original constraint tree missing
  taxa have been placed
}
\details{
  An analysis of the placement of each taxon not contained
  in the constraint tree is conducted. For each such taxon
  the edge lengths in the constraint tree are adjusted
  according to the parameter edge_scaling. If
  simple_edge_scaling is TRUE, the edges will have length 0
  if a taxon is never descendant from a tree and 1 if it is
  descendant in at least one tree. If FALSE, indicate the
  proportion of sampled trees in which the taxon is
  descendant from that edge.

  The output from this function is useful for checking that
  missing taxa are placed in appropriate positions relative
  to the original constraint tree.

  Note that this routine has not been optimised and slow
  (possibly unacceptable) performance is to be expected
  with large trees and/or large posterior samples.
}
\examples{
data(pastis_data_3_trees)

\dontrun{

# Check constraints for all missing taxa (takes ~20 seconds to run)
conch(pastis_data_3_trees[[1]], pastis_data_3_trees[[2]])

}

# Check constraints for missing taxon "a_4"
conch(pastis_data_3_trees[[1]], pastis_data_3_trees[[2]], species_set="a_4")
}

