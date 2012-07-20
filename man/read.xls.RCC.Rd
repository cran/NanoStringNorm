\name{read.xls.RCC}
\title{read.xls.RCC}
\alias{read.xls.RCC}
\description{A function to read the raw counts from the RCC excel spreadsheet output by the nCounter platform.}
\usage{
read.xls.RCC(xls, sheet = 1, perl, sample.id.row = 'File.Name')
}
\arguments{
  \item{xls}{The excel spreadsheet output by nCounter.  This should be a string with the path pointing to the desired file.}
  \item{sheet}{The worksheet that contains the raw counts.  Make sure to check that you are using the worksheet with the "raw" counts and not something that has been processed.  The name of the correct worksheet usually has "RCC" in it but not "norm".  The input is an integer and it defaults to 1 or the first worksheet. }
  \item{perl}{The path to a perl binary.  This does not need to be specified if perl is in your PATH}
\item{sample.id.row}{The row in the RCC file that contains the sample IDs to be used as column names for the count data.  Defaults to "File.Name" but in some cases "Sample.ID" or "Lane.ID" would be more appropriate.}
}

\value{
Returns a list with two components.  The first is the header information which contains sample IDs and diagnostic information on the quality of the samples.  The second is the count data and can be directly used in the input to NanoStringNorm.
}
\examples{
 
# directly import the nCounter output
path.to.xls.file <- system.file("extdata", "RCC_files", "RCCCollector1_rat_tcdd.xls", package = "NanoStringNorm");
NanoString.mRNA <- read.xls.RCC(x = path.to.xls.file, sheet = 1);

}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
