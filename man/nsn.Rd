\name{nsn}
\title{Methods for preprocessing NanoString nCounter data.}
\alias{nsn}
\description{This function can be used to normalize mRNA and miRNA expression data from the NanoString platform.  This is just a convenience wrapper of the NanoStrinNorm function}
\usage{
nsn(
	x, 
	...
	);
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an Excel spreadsheet.  If \emph{anno} is \emph{NA} then the first three columns must be c('Code.Class', 'Name', 'Accession') and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes. For convenience you can use the Excel import functions \emph{read.xls.RCC} to read directly from nCounter output files.  }
  \item{...}{Other NanoStringNorm arguments }
}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
