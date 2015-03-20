\name{read.markup.RCC}
\title{read.markup.RCC}
\alias{read.markup.RCC}
\description{A function to read and merge a set of sample specific RCC markup files.}
\usage{
read.markup.RCC(
 rcc.path = ".",
 rcc.pattern = "*.RCC|*.rcc",
 exclude = NULL, 
 include = NULL,
 nprobes = -1)
}
\arguments{
  \item{rcc.path}{The directory of the rcc files. Defaults to current working directory.}
  \item{rcc.pattern}{The file pattern used to search for the RCC files.  Defaults to the RCC extension}
  \item{exclude}{A list of files to ignore}
  \item{include}{A list of files to include}
  \item{nprobes}{(Optional) This defaults to the read.table default for nrows (-1). In some cases, a message follows the comment ("<") line and causes an error. In this case, the user can now define the number of rows (number of probes including controls and experimental) that should be read in.}
}

\value{
Returns a list with two components.  The first is the header information which contains sample IDs and diagnostic information on the quality of the samples.  The second is the count data and can be directly used in the input to NanoStringNorm.
}
\examples{

# read in all the rcc files in the current directory 
# data.raw <- read.markup.RCC();
# data.norm <- NanoStringNorm(data.raw);

}
\author{Daryl M. Waggott}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
