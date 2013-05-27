\name{create_constraint}
\alias{create_constraint}
\title{Constraint creator}
\usage{
  create_constraint(type, include, exclude = c(),
    name = NA, name_class = NA)
}
\arguments{
  \item{type}{constraint type either "hard","negative" or
  "partial"}

  \item{include}{the first list of taxa}

  \item{exclude}{the second list of taxa (not used for hard
  constraints)}

  \item{name}{the name for the constraint. If none
  specified they are automatically assigned}

  \item{name_class}{alternatively automatically incremented
  names}
}
\value{
  A string containing the constraint
}
\description{
  Creates mrBayes style constraints for the supplied list
  of species. Constraints are named according to either the
  supplied name or by the supplied name_class. If
  name_class is used each new constraint is incremented.
}

