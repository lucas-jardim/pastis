\name{collapse_cherry}
\alias{collapse_cherry}
\title{Collapse a cherry (may be a polytomy)}
\usage{
  collapse_cherry(tree, edge = NA, combine = NA)
}
\arguments{
  \item{tree}{An ape tree}

  \item{edge}{The edge to collapse (NA to choose the first
  cherry)}

  \item{combine}{For polytomies use this to choose which
  tips to combine. Specified as a list containing elements
  between 1 and the number of tips in the polytomy. Tips
  are numbered according to their order in tree$edge. Omit
  or set to NA to combine all tips.}
}
\value{
  A list containing: 1) the tree with a collapsed cherry 2)
  an index to the edge that was collapsed
}
\description{
  Collapse a cherry (may be a polytomy)
}

