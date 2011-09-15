\name{get.trait.summary.stats}
\title{get.trait.summary.stats}
\alias{get.trait.summary.stats}
\description{A function to get summary information on genes in a NanoString dataset. }
\usage{
get.trait.summary.stats(x, anno, log, traits)
}
\arguments{
  \item{x}{A NanoString RCC worksheet.  If anno is not specified then the first three columns must be c("Code.Class","Name","Accession").  The remaining columns are samples and rows are genes. }
  \item{anno}{An alternative and flexible annotation dataset.  The columns c("Code.Class","Name","Accession") do not need to be included in "x" if anno is specified.}
  \item{log}{Has the data been logged (log2).  Raw NanoString data is not in log space.  This will only effect the calculation of Fold-change.}
  \item{traits}{A vector or matrix of binary traits such as normal/tumour status.  Each trait may only contain 1,2 or NA similar to the numeric status of factors.  T-tests p-values and Fold-Change in terms of the "2" category are presented.  The results can be displayed using build in plotting functions.}
By default the function returns a list of objects including pre-processed counts, specified model, sample and gene summaries.  The sample summary holds the calculated normalization factors for evaluation of problem samples.  The gene summary includes means, coefficient of variation and differential expression statistics for any traits. 
}
\value{
A matrix including t-test p-values and Fold-Changes for all genes and traits.
}
\author{Daryl M. Waggott}

\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
