\name{get.batch.effects}
\title{get.batch.effects}
\alias{get.batch.effects}
\description{A function to test if there is an association between binary traits and the sample normalization factors.  Traits can include different run times, sample conditions, plate etc.  This function is called internally from the get.trait.summary.stats function. }
\usage{
get.batch.effects(x, anno, log, traits, sample.summary.stats)
}
\arguments{
  \item{x}{A NanoString RCC worksheet.  If anno is not specified then the first three columns must be c("Code.Class","Name","Accession").  The remaining columns are samples and rows are genes. }
  \item{anno}{An alternative and flexible annotation dataset.  The columns c("Code.Class","Name","Accession") do not need to be included in "x" if anno is specified.}
  \item{log}{Has the data been logged (log2).  Raw NanoString data is not in log space.  This will only effect the calculation of Fold-change.}
  \item{traits}{A vector or matrix of binary traits such as normal/tumour status.  Each trait may only contain 1,2 or NA similar to the numeric status of factors.  T-tests p-values and Fold-Change in terms of the "2" category are presented.  The results can be displayed using build in plotting functions.}
  \item{sample.summary.stats}{a data.frame of sample summary stats such as mean, sd, background level and normalization factors.}
}
\value{
A matrix where each row has results for a one summary statistic vs one trait.
}
\author{Daryl M. Waggott}

\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
