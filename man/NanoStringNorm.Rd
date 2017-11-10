\name{NanoStringNorm}
\title{Methods for preprocessing NanoString nCounter data}
\alias{NanoStringNorm}
\description{This function can be used to normalize mRNA and miRNA expression data from the NanoString platform.}
\usage{
NanoStringNorm(
	x, 
	anno = NA,
	header = NA,
	Probe.Correction.Factor = 'adjust',
	CodeCount = 'none',
	Background = 'none',
	SampleContent = 'none',
	OtherNorm = 'none',
	CodeCount.summary.target = NA,
	SampleContent.summary.target = NA,
	round.values = FALSE,
	is.log = FALSE,
	take.log = FALSE,
	return.matrix.of.endogenous.probes = FALSE,
	traits = NA,
	predict.conc = FALSE,
	verbose = TRUE,
	genes.to.fit,
	genes.to.predict,
	guess.cartridge = TRUE,
	...
	);
}
\arguments{
  \item{x}{The data used for Normalization.  This is typically the raw expression data as exported from an Excel spreadsheet.  If \emph{anno} is \emph{NA} then the first three columns must be c('Code.Class', 'Name', 'Accession') and the remaining columns refer to the samples being analyzed.  The rows should include all control and endogenous genes. For convenience you can use the Excel import functions \emph{read.xls.RCC} to read directly from nCounter output files.  }
  \item{anno}{Alternatively, \emph{anno} can be used to specify the first three annotation columns of the expression data.  If anno used then it is assumed that \emph{x} does not contain these data.  Anno allows flexible inclusion of alternative annotation data.  The only requirement is that it includes the \emph{Code.Class}  and \emph{Name} columns. Code.Class refers to gene classification i.e. Positive, Negative, Housekeeping or Endogenous gene.  For some Code Sets you will need to manually specify Housekeeping genes.}
  \item{header}{The sample meta data at the top of the RCC excel worksheet.  The data typically spans rows 4-18 and there should be the same number of columns as samples.  This is imported automaically when you use the \emph{read.xls.RCC}.  Values in the header will be parsed for normalization diagnostics i.e. is expression associated with a specific cartridge.}
  \item{Probe.Correction.Factor}{An adjustment factor to be applied at the probe level prior to any normalization.  Specify 'filter' if you would like to remove any flagged genes.  The default to 'adjust' which uses the probe correction factors concatenated to the gene names.  See details.}
  \item{CodeCount}{The method used to normalize for technical assay variation.  The options are \emph{none} and \emph{sum} and \emph{geo.mean} (geometric mean).  The method adjusts each sample based on its relative value to all samples.  The geometric mean may be less susceptible to extreme values.  The CodeCount normalization is applied first and is considered the most fundamental Normalization step for technical variation.}
  \item{Background}{The method used to estimate the background count level.  Background is calculated based on negative controls and the options are \emph{none, mean, mean.2sd, max}.  Mean is the least is the least conservative, while max and mean.2sd are the most robust to false positives.  The calculated background is subtracted from each sample.  Background is calculated after code count adjustment.}
  \item{SampleContent}{The method used to normalize for sample or RNA content i.e. pipetting fluctuations.  The options are \emph{none, housekeeping.sum, housekeeping.geo.mean, total.sum, low.cv.geo.mean, top.mean} and \emph{top.geo.mean}.  The Housekeeping options require a set of annotated genes.  For RNA datasets the genes can be specified by editing the Code.Class field to equal 'Housekeeping' or 'Control'.  miRNA Code Sets generally already have annotated Housekeeping genes.  Note that the Housekeeping genes on miRNA Code Sets are messenger RNAs and therefore undergo a different laboratory process, specifically no ligation.  If using Housekeeping genes then 'housekeeping.geo.mean' is recommended.  Alternatively, you can use the method \emph{total.sum} which uses the sum of the all the genes or \emph{top.geo.mean} which uses the geometric mean of the top 75 expressed.  \emph{low.cv.geo.mean} has recently been added to reduce the influence of outliers and normalizing to trait associated genes.  This method simply selects genes that expressed with low coefficients of variation.  Sample Content adjustment is applied after Code Count and Background correction.  }
  \item{OtherNorm}{Some additional normalization methods.  Options include 'none', 'quantile' and 'zscore'.  OtherNorm is applied after CodeCount, Background, SampleContent normalizations.  This is probably the largest source of inter sample variation in most studies.}
  \item{CodeCount.summary.target}{The expected summary value for CodeCount positive control probes. By default this is calculated as the mean of summary values across all samples. This target value is then used to determine the normalization factor of each sample when a CodeCount method is applied. User-assignment of the target value produces sample-independent CodeCount normalization, allowing for new samples to be added without changing the normalization factors of existing samples.}
  \item{SampleContent.summary.target}{The expected summary value for SampleContent housekeeping control probes. By default this is calculated as the mean of summary values across all samples. This target value is then used to determine the normalization factor of each sample when a SampleContent method is applied. User-assignment of the target value produces sample-independent SampleContent normalization, allowing for new samples to be added without changing the normalization factors of existing samples. Note that by definition 'low.cv.geo.mean', 'top.mean', and 'top.geo.mean' methods are not compatible with sample-independent normalization.}
  \item{round.values}{Should the values be rounded to the nearest absolute count or integer.  This simplifies interpretation if taking the log by removing values 0-1 which result in negatives.}
  \item{is.log}{Is the data already in logspace.  This is recommended if applying NSN to PCR type microarray data.  In this scenario the geometric log is altered and negative values are not flagged.}
  \item{take.log}{Should a log2 transformation of the data be used.  Taking the log may help distributional assumptions and consistency with PCR based methods for calculating fold change.  Note, this is not a flag on if the log has already been taken.}
  \item{return.matrix.of.endogenous.probes}{If true a matrix of normalized Endogenous code counts is returned.  Annotation columns and control probes are removed.  This can be useful if you the output is being used directly for downstream analysis.  By default it is FALSE and list of objects including descriptive and diagnostics are returned.}
  \item{traits}{A vector or matrix of phenotypes, design or batch effect variables.  For example tumour status, trx regimen, age at FFPE fixation, RNA quality scores and sample plate.  Note, that nCounter design covariates such as 'cartridge' and 'date' can be collected from the header of the RCC Excel spreadsheet.  At this time each trait should be binary and may only contain 1,2 or NA similar to the numeric coding of factors.  T-tests p-values and Fold-Change are presented in terms of the '2' category.  The results can be displayed using built in plotting functions.}
  \item{predict.conc}{Should the predicted concentration values be returned.  Defaults to FALSE.  The predicted concentrations are based on a fitted model between the observed counts and the expected concentration (fM).  Counts are replaced by predicted concentrations from a model using each samples normalized data. In addition a final column is added which predicts the concentration using the mean of all samples.}
  \item{verbose}{Output comments on run-status}
  \item{genes.to.fit}{The set of genes used to generate the normalization parameters.  You can specify a vector of code classes or gene names as indicated in the annotation.  Alternatively, you can use descriptors such as 'all', 'controls', 'endogenous'.  For most methods the model will be fit and applied to the endogenous genes.  In some cases such as vsn you may want to vary this approach.}
  \item{genes.to.predict}{The set of genes that the parameters or model is applied to.}
  \item{guess.cartridge}{Should the cartridge be estimated based on the order of the samples.}
  \item{...}{The ellipses are included to allow flexible use of parameters that are required by OtherNorm methods.  For example the use of \emph{strata or calib} parameters documented in the vsn package.  }
}
\details{
The code is based on the NanoString analysis guidelines (see references).  The function allows normalization of both mRNA and miRNA NanoString expression data.  The order of the methods is fixed but the use of the method \emph{none} and multiple iterations of the functions allows flexibility.

Note.  Poorly assayed samples could negatively influence the normalization of the remaining data.  Prior to normalization check that the binding density is not less than .04 and the number of total counts/FOV is not much less than 1500. 

In the newer RCC format (Spring 2012) the probe correction factors are concatenated to the gene name via a pipe symbol.  By default these are parsed and corrected for.  In the older RCC format the 'Name' column of the RCC worksheet sometimes flags certain probes with the message '(+++ See Message Below)'.  If this is the case a 'Readme' document including probe level adjustment factors should have been supplied by your Microrray center.  This file must be edited into a tabular file and specified in the Probe.Correction.Factor argument.  The function will fail with an error if warnings are detected and no probe levels correction factor is supplied. Upon correction any warning messages will be stripped from the raw data file.  This background is probe-specific and all background subtraction should be completed prior to normalization.  The number of counts to be subtracted for a given probe is determined by multiplying a correction factor by the number of counts observed for the 128fM positive control in the same lane.

}
\value{
By default the function returns a list of objects including pre-processed counts, specified model, sample and gene summaries.  The sample summary holds the calculated normalization factors for evaluation of problem samples.  The gene summary includes means, coefficient of variation and differential expression statistics for any traits. 
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

### 1 Normalize mRNA and output a matrix of normalized counts #################

# load the NanoString.mRNA dataset
data(NanoString);

# alternatively you can import the dataset
# path.to.xls.file <- system.file("extdata", "RCC_files", "RCCCollector1_rat_tcdd.xls", 
# 	package = "NanoStringNorm");
# NanoString.mRNA <- read.xls.RCC(x = path.to.xls.file, sheet = 1);

# specifiy housekeeping genes in annotation 
NanoString.mRNA[NanoString.mRNA$Name \%in\% 
	c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';

# normalize
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	anno = NA,
	CodeCount = 'geo.mean',
	Background = 'mean',
	SampleContent = 'housekeeping.geo.mean',
	round.values = TRUE,
	take.log = TRUE,
	return.matrix.of.endogenous.probes = TRUE
	);

### 2 include a trait for differential expression and batch effect evaluation  #
### 	A list of diagnostic objects is output.                                #

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

# Split the input into separate annotation and count input  ###################

NanoString.mRNA.anno <- NanoString.mRNA[,c(1:3)];
NanoString.mRNA.data <- NanoString.mRNA[,-c(1:3)];

# Normalize
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA.data,
	anno = NanoString.mRNA.anno,
	CodeCount = 'sum',
	Background = 'mean.2sd',
	SampleContent = 'top.geo.mean',
	round.values = TRUE,
	take.log = TRUE,
	traits = trait.strain,
	return.matrix.of.endogenous.probes = FALSE
	);

### 3 plot results ############################################################

# Plot all the plots as PDF report.  See help on Plot.NanoStringNorm for examples
pdf('NanoStringNorm_Example_Plots_All.pdf');
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = 'all'
	);
dev.off();

}
\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
