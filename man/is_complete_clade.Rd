\name{is_complete_clade}
\alias{is_complete_clade}
\title{Clade completeness check}
\usage{
  is_complete_clade(set, clades)
}
\arguments{
  \item{set}{A vector of taxa}

  \item{clades}{A list of clades, each element contains a
  vector of taxa belonging to that clade}
}
\value{
  logical, whethe the clades were completely represented
}
\description{
  Determines all clades represented by a set of taxa and
  checks whether all species in those clades are in that
  set of taxa
}

