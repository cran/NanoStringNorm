\name{get.sample.summary.stats}
\title{get.sample.summary.stats}
\alias{get.sample.summary.stats}
\description{A function to get basic sample summary information including mean, sd and coefficient of variation. }
\usage{
get.sample.summary.stats(x, anno)
}
\arguments{
  \item{x}{A NanoString RCC worksheet.  If anno is not specified then the first three columns must be c("Code.Class","Name","Accession").  The remaining columns are samples and rows are genes. }
  \item{anno}{An alternative and flexible annotation dataset.  The columns c("Code.Class","Name","Accession") do not need to be included in "x" if anno is specified.}
}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
