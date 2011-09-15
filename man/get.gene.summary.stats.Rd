\name{get.gene.summary.stats}
\title{get.gene.summary.stats}
\alias{get.gene.summary.stats}
\description{A function to get basic gene summaries such as mean, sd, and coefficient of variation}
\usage{
get.gene.summary.stats(x, anno = NA)
}
\arguments{
  \item{x}{A NanoString RCC worksheet.  If anno is not specified then the first three columns must be c("Code.Class","Name","Accession").  The remaining columns are samples and rows are genes. }
  \item{anno}{An alternative and flexible annotation dataset.  The columns c("Code.Class","Name","Accession") do not need to be included in "x" if anno is specified.}
}
\value{
A data.frame with mean, standard deviation and coefficient of variation for each gene on the log scale.
}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
