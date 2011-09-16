\name{predict.concentration}
\title{predict.concentration}
\alias{predict.concentration}
\description{A function to predict the concentration of endogenous genes from their counts.  The model is based on the relationship between the positive control counts and their known concentration.}
\usage{
predict.concentration(x, anno, log, verbose = TRUE);
}
\arguments{
  \item{x}{A NanoString RCC worksheet.  If anno is not specified then the first three columns must be c("Code.Class","Name","Accession").  The remaining columns are samples and rows are genes. }
  \item{anno}{An alternative and flexible annotation dataset.  The columns c("Code.Class","Name","Accession") do not need to be included in "x" if anno is specified.}
  \item{log}{Has the data been logged (log2).  Raw NanoString data is not in log space.  This will only effect the calculation of Fold-change.}
By default the function returns a list of objects including pre-processed counts, specified model, sample and gene summaries.  The sample summary holds the calculated normalization factors for evaluation of problem samples.  The gene summary includes means, coefficient of variation and differential expression statistics for any traits. 
  \item{verbose}{Output run-time status messages}
}
\value{
A matrix including t-test p-values and Fold-Changes for all genes and traits.
}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
