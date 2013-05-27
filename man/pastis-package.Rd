\docType{package}
\name{pastis-package}
\alias{pastis}
\alias{pastis-package}
\title{Phylogenetic Assembly with Soft Taxonomic Inferences.}
\description{
  Phylogenetic Assembly with Soft Taxonomic Inferences.
}
\details{
  \tabular{ll}{ Package: \tab pastis\cr Type: \tab
  Package\cr Version: \tab 0.1\cr Date: \tab 2013-05-27\cr
  Depends: \tab caper (>= 0.5.0), ape (>= 3.0.0)\cr
  Encoding: \tab UTF-8\cr License: \tab GPL (>= 3)\cr
  LazyLoad: \tab yes\cr URL: \tab
 \cr }

  A pre-processor for mrBayes that assimilates sequences,
  taxonomic information and tree constraints as per xxx.
  The main functions of interest for most users will be
  pastis_simple, pastis_main and conch. The main analysis
  is conducted with pastis_simple or pastis_main followed
  by a manual execution of mrBayes (>3.2). The placement of
  taxa not contained in the tree constraint can be
  investigated using conch.
}
\author{
  Klaas Hartmann \email{klaas.hartmann@utas.edu.au}, 
  Gavin Thomas \email{gavin.thomas@sheffield.ac.uk}
}

\keyword{package}

