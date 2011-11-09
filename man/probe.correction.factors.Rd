\name{probe.correction.factor.normalization}
\title{probe.correction.factor.normalization}
\alias{probe.correction.factor.normalization}
\description{This function can be used to normalize mRNA and miRNA expression data from the NanoString platform.}
\usage{
probe.correction.factor.normalization(x, anno, Probe.Correction.Factor, verbose = TRUE);
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an XLS spreadsheet.  If anno is NA then the first three columns must be c('Code.Class', 'Name', 'Accession') and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes.}
  \item{anno}{Alternatively anno can be used to specify the first three annotation columns of the expression data.  If anno used then it assumed that 'x' does not contain these data.  Anno allows flexible inclusion of alternative annotation data.  The only requirement is that it includes the "Code.Class" and "Name" columns. Code.Class refers to gene classification i.e. Positive, Negative, Housekeeping or Endogenous gene.}
  \item{Probe.Correction.Factor}{An adjustment factor to be applied at the probe level prior to any normalization.  The data should be a matrix with one 2 columns, one for gene Name and another for probe CorrectionFactor).  Specify "filter" if you would like to remove any flagged genes.  See details.}
  \item{verbose}{Output run-time status messages}
}
\details{
Poorly assayed samples could negatively influence the normalization of the remaining data.  Prior to normalization check that the binding density is not less than .05 and the number of total counts/FOV is not much less than 1500. 
The "Name" column of the RCC worksheet sometimes flags certain probes with the message "(+++ See Message Below)".  If this is the case a "README" document including probe level adjustment factors should have been supplied by your Microrray center.  This file must be edited into a tabular file and specified as the argument to the Probe.Correction.Factor parameter.  The function will fail with an error if warnings are detected and no probe levele correction factor is supplied. Upon correction any warning messages will be stripped from the raw data file.  This background is probe-specific and all background subtraction should be completed prior to normalization.  The number of counts to be subtracted for a given probe is determined by multiplying a correction factor by the number of counts observed for the 128 fM positive control in the same lane.
}
\references{
See NanoString website for PDFs on analysis guidelines:
	\url{http://www.nanostring.com/support/prod-lit/}

The NanoString assay is described in the paper:
	Fortina, P. & Surrey, S. Digital mRNA profiling. Nature biotechnology 26, 293-294 (2008).
}
\author{Daryl M. Waggott}
\note{Future updates to include more informative diagnostic plots, trait specific normalization and replicate merging.}
\examples{
# INSERT EXAMPLES HERE
}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
