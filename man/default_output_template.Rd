\name{default_output_template}
\alias{default_output_template}
\title{A sample output template for pastis_main and pastis_simple}
\usage{
  default_output_template()
}
\value{
  The default output template
}
\description{
  This is the default output template filled in by
  pastis_main and pastis_simple to create the input file
  for mrBayes.
}
\details{
  If you want to change the parameters in the mrBayes nexus
  file created by pastis then it is easier to create your
  own template rather than editing the created nexus files
  manually (which you would have to do again if you reran
  pastis).

  This is the default template used by pastis and is a good
  starting point for creating your own template. To have a
  look at this template try:

  \verb{ cat(default_output_template()) }

  Once you understand the format, you can create your own
  template as a string and pass it to pastis.
}

