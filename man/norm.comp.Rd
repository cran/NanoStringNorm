\name{norm.comp}
\title{norm.comp}
\alias{norm.comp}
\description{Compare normalizaton methods using signal to noise CV and replicates ICC.}
\usage{
norm.comp(x, anno, replicates = NULL, CodeCount.methods = c('none', 'sum', 'geo.mean'), Background.methods = c('none','mean', 'mean.2sd','max'), 
SampleContent.methods = c('none','housekeeping.sum', 'housekeeping.geo.mean', 'total.sum','top.mean', 'top.geo.mean', 'low.cv.geo.mean'),
OtherNorm.methods = c('none','quantile','zscore', 'rank.normal', 'vsn'),  histogram = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an Excel spreadsheet.  If anno is NA then the first three columns must be \emph{c('Code.Class', 'Name', 'Accession')} and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes.}
  \item{anno}{Alternatively, anno can be used to specify the first three annotation columns of the expression data.  If anno used then it assumed that 'x' does not contain these data.  Anno allows flexible inclusion of alternative annotation data.  The only requirement is that it includes the 'Code.Class' and 'Name' columns. Code.Class refers to gene classification i.e. Positive, Negative, Housekeeping or Endogenous gene.}
  \item{replicates}{A vector of IDs indicating what samples are replicates.  Replicate samples need to have the same ID.}
  \item{CodeCount.methods}{vector of methods to compare}
  \item{Background.methods}{vector of methods to compare}
  \item{SampleContent.methods}{vector of methods to compare}
  \item{OtherNorm.methods}{vector of methods to compare}
  \item{histogram}{logical if histogram is to be output}
  \item{verbose}{logical if logging should be output}
}
\details{
A data.frame is returned with a single row for each combination of methods.  For each method the coefficient of variation is reported for positive controls, housekeeping and endogenous genes.  If replicates are included then intra-class correlation (ICC) values are reported.  One would expect that minimizing variation in controls and replicates is a reasonable estimate of normalization performance.  Ideally, this should be performed on a pilot study.  Future extensions will include comparing methods via cross validation on trait differences in order to maximize information content.
}
\author{Daryl M. Waggott}

\examples{

# load the NanoString.mRNA dataset
data(NanoString);


# specifiy housekeeping genes in annotation
NanoString.mRNA[NanoString.mRNA$Name \%in\% c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';

# strain x experimental condition i.e. replicate.
# this is only a small subset of the original data used for the plot
biological.replicates <- c("HW_1.5_0","HW_1.5_0","HW_1.5_0","HW_1.5_100","HW_1.5_100","HW_1.5_100","HW_6_100","HW_6_100","HW_3_100",
"HW_3_100","HW_3_100","HW_3_100","LE_19_0","LE_19_0","LE_19_0","LE_96_0","LE_96_0","LE_96_0","HW_10_100","HW_10_100","HW_10_100",
"HW_10_100","HW_6_100","HW_6_100","HW_96_0");

norm.comp.results <- norm.comp(
x = NanoString.mRNA,
replicates = biological.replicates,
CodeCount.methods = 'none',
Background.methods = 'none',
SampleContent.methods = c('none','housekeeping.sum', 'housekeeping.geo.mean','top.mean', 'top.geo.mean'),
OtherNorm.methods = 'none',
verbose = FALSE
);

}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
