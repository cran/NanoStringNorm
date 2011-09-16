\name{NanoStringNorm}
\title{NanoStringNorm}
\alias{NanoStringNorm}
\description{This function can be used to normalize mRNA and miRNA expression data from the NanoString platform.}
\usage{
NanoStringNorm(x, anno = NA, Probe.Correction.Factor = 'none', CodeCount = 'none', Background = 'none', SampleContent = 'none', otherNorm = 'none', round.values = FALSE, log = FALSE, return.matrix.of.endogenous.probes = FALSE, traits = NA, verbose = TRUE);
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an XLS spreadsheet.  If anno is NA then the first three columns must be c('Code.Class', 'Name', 'Accession') and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes.}
  \item{anno}{Alternatively anno can be used to specify the first three annotation columns of the expression data.  If anno used then it assumed that 'x' does not contain these data.  Anno allows flexible inclusion of alternative annotation data.  The only requirement is that it includes the 'Code.Class' and 'Name' columns. Code.Class refers to gene classification i.e. Positive, Negative, Housekeeping or Endogenous gene.}
  \item{Probe.Correction.Factor}{An adjustment factor to be applied at the probe level prior to any normalization.  Specify 'filter' if you would like to remove any flagged genes.  See details.}
  \item{CodeCount}{The method used to normalize for the majority of technical assay variation.  The options are 'none' and 'sum' and 'geo.mean' (geometric mean).  The method adjusts each sample based on its relative value to all samples.  The geometric mean may be less susceptible to extreme values.  The CodeCount normalization is applied first and is considered the most fundamental Normalization step.}
  \item{Background}{The method used to estimate the background level.  Background is calculated based on negative controls and the options are 'none', 'mean', 'mean.2sd', 'max'.  Mean is the least and max is the most conservative and robust to false positives.  The calculated background is subtracted from each sample.  Background is calculated after code count adjustment.}
  \item{SampleContent}{The method used to normalize for sample or RNA content i.e. pipetting fluctuations.  The options are 'none', 'housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'top.mean' and 'top.geo.mean'.  The housekeeping options require a set of annotated genes.  For RNA datasets the genes can be specified by editing the Code.Class field to equal 'Housekeeping' or 'Control' and the recommended method is 'housekeeping.geo.mean'.  For miRNA datasets the housekeeping genes are already annotated and the 'housekeeping.mean' method is recommended.  'total' uses the sum of the all the genes while 'top' uses the mean of the top 100 expressed genes which minimizes the impact of undetected genes.  The option 'top.geo.mean' is an alternative method which ranks RNAs based on the sum of all samples and then takes the geometric mean of the top 75.  Sample Content adjustment is applied after Code Count and Background correction.}
  \item{otherNorm}{Some additional normalization methods.  Options include 'none', 'quantile' and 'zscore'.  otherNorm is applied after CodeCount, Background, SampleContent normalizations.}
  \item{round.values}{Should the values be rounded to the nearest absolute count or integer.  This simplifies interpretation if taking the log by removing values 0-1 which result in negatives.}
  \item{log}{Should a log2 transformation of the data be used.  Taking the log may help distributional assumptions and consistency with PCR based methods of calculating fold change.}
  \item{return.matrix.of.endogenous.probes}{If true a matrix of normalized endogenous code counts is returned.  Annotation columns and control probes are removed.  This could be useful if you the output is being used directly for downstream analysis.  By default this is FALSE and list of objects including descriptive and diagnostics.}
  \item{traits}{A vector or matrix of binary traits such as normal/tumour status.  Each trait may only contain 1,2 or NA similar to the numeric status of factors.  T-tests p-values and Fold-Change in terms of the '2' category are presented.  The results can be displayed using build in plotting functions.}
  \item{verbose}{Output comments on run-status}
}
\details{
The code is based on the NanoString analysis guidelines(see references).  The function allows normalization of both mRNA and miRNA NanoString expression data.  The order of the methods is fixed but the use of NA and multiple iterations of the functions allows flexibility.

Note.  Poorly assayed samples could negatively influence the normalization of the remaining data.  Prior to normalization check that the binding density is not less than .05 and the number of total counts/FOV is not much less than 1500. 

The 'Name' column of the RCC worksheet sometimes flags certain probes with the message '(+++ See Message Below)'.  If this is the case a 'Readme' document including probe level adjustment factors should have been supplied by your Microrray center.  This file must be edited into a tabular file and specified in the Probe.Background.Correction argument.Factor.  The function will fail with an error if warnings are detected and no probe levele correction factor is supplied. Upon correction any warning messages will be stripped from the raw data file.  This background is probe-specific and all background subtraction should be completed prior to normalization.  The number of counts to be subtracted for a given probe is determined by multiplying a correction factor by the number of counts observed for the 128fM positive control in the same lane.

}
\value{
By default the function returns a list of objects including pre-processed counts, specified model, sample and gene summaries.  The sample summary holds the calculated normalization factors for evaluation of problem samples.  The gene summary includes means, coefficient of variation and differential expression statistics for any traits. 
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

### 1 Normalize mRNA and output a matrix of normalized counts ###

# load data
data(NanoString);

# specifiy housekeeping genes in annotation 
NanoString.mRNA[NanoString.mRNA$Name \%in\% c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';

# normalize
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	anno = NA,
	CodeCount = 'geo.mean',
	Background = 'mean',
	SampleContent = 'housekeeping.geo.mean',
	round.values = TRUE,
	log = TRUE,
	return.matrix.of.endogenous.probes = TRUE
	);

### 2 include a trait for differential expression and batch effect evaluation.  A list of diagnostic objects is output. ###

# Define a traits based on strain
sample.names <- names(NanoString.mRNA)[-c(1:3)];
strain1 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain1[grepl('HW',sample.names)] <- 2;
strain2 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain2[grepl('WW',sample.names)] <- 2;
strain3 <- rep(1, times = (ncol(NanoString.mRNA)-3));
strain3[grepl('LE',sample.names)] <- 2;
trait.strain <- data.frame(
	row.names = sample.names,
	strain1 = strain1,
	strain2 = strain2,
	strain3 = strain3
	);

# Normalize
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	anno = NA,
	CodeCount = 'sum',
	Background = 'mean.2sd',
	SampleContent = 'top.geo.mean',
	round.values = TRUE,
	log = TRUE,
	traits = trait.strain,
	return.matrix.of.endogenous.probes = FALSE
	);

### 3 plot results ###

# Plot all the plots as PDF report.  See help on Plot.NanoStringNorm for further examples.
pdf('NanoStringNorm_Example_Plots_All.pdf')
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = 'all'
	);
dev.off()

}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
