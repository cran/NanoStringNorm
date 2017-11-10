\name{other.normalization}
\title{other.normalization}
\alias{other.normalization}
\description{This function can be used to normalize mRNA and miRNA expression data from the NanoString platform.}
\usage{
other.normalization(
 x,
 anno,
 OtherNorm = 'none',
 verbose = TRUE,
 genes.to.fit,
 genes.to.predict,
 ...);
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an Excel spreadsheet.  If anno is NA then the first three columns must be \emph{c('Code.Class', 'Name', 'Accession')} and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes.}
  \item{anno}{Alternatively, anno can be used to specify the first three annotation columns of the expression data.  If anno used then it assumed that 'x' does not contain these data.  Anno allows flexible inclusion of alternative annotation data.  The only requirement is that it includes the 'Code.Class' and 'Name' columns. Code.Class refers to gene classification i.e. Positive, Negative, Housekeeping or Endogenous gene.}
  \item{OtherNorm}{Some additional normalization methods.  Options \strong{currently} include \emph{'none', 'zscore', 'rank.normal'} and \emph{'vsn'}.  OtherNorm is applied after CodeCount, Background, SampleContent normalizations, therefore you may want to set them to 'none'.}
  \item{verbose}{Enable run-time status messages}
  \item{genes.to.fit}{The set of genes used to generate the normalization parameters.  You can specify a vector of code classes or gene names as indicated in the annotation.  Alternatively, you can use descriptors such as 'all', 'controls', 'endogenous'.  For most methods the model will be fit and applied to the endogenous genes.  vsn defaults to 'all' genes.}
  \item{genes.to.predict}{The set of genes that the parameters or model is applied to.}
  \item{...}{The ellipses are included to allow flexible use of parameters that are required by sub normalization methods.  For example the use of \emph{'strata'} or \emph{'calib'} parameters documented in the vsn package.  }
}
\details{
The methods used for OtherNorm are designed to be extensible to alternate and evolving NanoString pre-processing analysis.  These can be combined with standard  CodeCount, Background, SampleContent methods (i.e. positive, negative and housekeeping controls).  
\cr\cr
zscore is simply scaling each sample to have a mean 0 and standard deviation of 1.  Similarly, rank.normal uses the quantiles or ranks to transform to the normal distribution with mean 0 and standard deviation of 1.  Both methods are helpful when comparing multiple datasets or platforms during meta or joint analysis due to the abstraction of effects sizes or expression scale.
\cr\cr
quantile normalization is based on pegging each sample to an empirically derived distribution.  The distribution is simply the mean expression across all samples for each gene rank.
\cr\cr
vsn normalization is well documented and additional parameterization can be read about via it's package documentation.  One interesting feature is that the affine transformation can be turned off leaving just the generalized log2 (glog2) transformation \emph{calib = 'none'}.  One can thereby use the existing NanoString methods to shift and scale the count data and then apply glog2 transformation to help stabilize the variance.
}
\references{
See NanoString website for PDFs on analysis guidelines:
	\url{https://www.nanostring.com/support/product-support/support-documentation}

The NanoString assay is described in the paper:
	Fortina, P. & Surrey, S. Digital mRNA profiling. Nature biotechnology 26, 293-294 (2008).
}
\author{Daryl M. Waggott}
\note{Future updates to include more informative diagnostic plots, trait specific normalization and replicate merging.}
\examples{

# load the NanoString.mRNA dataset
data(NanoString);

# specify housekeeping genes in annotation 
NanoString.mRNA[NanoString.mRNA$Name \%in\% 
	c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';

# z-value transformation. scale each sample to have a mean 0 and sd
# by default all the other normalization methods are 'none'
# you cannot apply a log because there are negative values
# good for meta-analysis and cross platform comparison abstraction of effect size
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	OtherNorm = 'zscore',
	return.matrix.of.endogenous.probes = TRUE
	);


# inverse normal transformation. use quantiles to transform each sample to 
# the normal distribution
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	OtherNorm = 'rank.normal',
	return.matrix.of.endogenous.probes = TRUE
	);

# quantile normalization.  create an empirical distribution based on the 
# median gene counts at the same rank across sample.  then transform each
# sample to the empirical distribution.
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	OtherNorm = 'quantile',
	return.matrix.of.endogenous.probes = FALSE
	);

if(require(vsn)) { # only run if vsn installed
# vsn.  apply a variance stabilizing normalization.
# fit and predict the model using 'all' genes 
#	i.e. 'controls' and 'endogenous', this is the default
# note this is just a wrapper for the vsn package.
# you could even add strata for the controls vs the 
# 	endogenous to review systematic differences

NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	OtherNorm = 'vsn',
	return.matrix.of.endogenous.probes = FALSE,
	genes.to.fit = 'all',
	genes.to.predict = 'all'
	);

# vsn. this time generate the parameters (fit the model) on the 'controls' and apply 
#	(predict) on the endogenous
# alternatively you may want to use the endogenous probes for both fitting and predicting.
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	OtherNorm = 'vsn',
	return.matrix.of.endogenous.probes = FALSE,
	genes.to.fit = 'controls',
	genes.to.predict = 'endogenous'
	);

# vsn. apply standard NanoString normalization strategies as an alternative to the vsn 
#	affine transformation.
# this effectively applies the glog2 variance stabilizing transformation to the 
#	adjusted counts
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	CodeCount = 'sum',
	Background = 'mean',
	SampleContent = 'top.geo.mean',
	OtherNorm = 'vsn',
	return.matrix.of.endogenous.probes = FALSE,
	genes.to.fit = 'endogenous',
	genes.to.predict = 'endogenous',
	calib = 'none'
	);
} # end vsn

}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
