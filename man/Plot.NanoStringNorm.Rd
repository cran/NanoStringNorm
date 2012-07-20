\name{Plot.NanoStringNorm}
\title{Plot.NanoStringNorm}
\alias{Plot.NanoStringNorm}
\description{Do some diagnostic and descriptive plots}
\usage{
Plot.NanoStringNorm(
	x,
	plot.type = 'norm.factors',
	samples = NA,
	genes = NA,
	label.best.guess = TRUE,
	label.ids = NA,
	title = TRUE,
	col = NA
	);

}
\arguments{
  \item{x}{An object of class NanoStringNorm i.e. the list not matrix output of NanoStringNorm}
  \item{plot.type}{A vector of the different plots that you are interested in.  Currently the following are implemented: normalization factors 'norm.factors', coefficient of variance 'cv' and mean vs sd 'mean.sd', 'volcano', 'missing', 'RNA.estimates', 'batch.effects' and 'positive.controls'.  To plot all available plots you can use 'all'. }
  \item{samples}{List of samples to use.  Not yet implemented.}
  \item{genes}{List of genes to use.  Not yet implemented.}
  \item{label.best.guess}{Labels will automatically be added for samples and genes that are either outliers (~3sd from mean) or the most significant.  The default is TRUE.}
  \item{label.ids}{Explicitly label points.  To use this feature supply a list with with two elements named genes and samples.  By default this is set to NA.}
  \item{title}{Should the plot titles be included.  By default they are included but this can be turned off if generating plots for publication or presentation.}
  \item{col}{Alternate colours for plotting.  Input should be a vector of two rgb values.  The default colours are orange and green.}
}
\value{
The output is single or series of plots.  Plots focus on either genes or samples.  Due to variation in sample quality it is a good idea to review data in order to determine outliers.  Outliers could effect normalization of other samples and bias results of downstream statistical analysis.  Systematic batch effects in terms of cartridge, sample type or covariates should also be evaluated to minimize confounding results.
\cr \cr
The mean vs standard deviation plot for each gene is a good indication of the noisiness of controls genes.  Also, candidate housekeeping genes can be identified by having a high mean and very low variance.
\cr \cr
The coefficient of variation density plot shows the change in signal to noise pre and post normalization.  
\cr\cr
The volcano plots show the relationship between p-value and fold-change for traits supplied in the normalization step.  These can be phenotype's of interest or study design based for evaluating batch effects.  Orange dots have a fold-change greater than 2 and are sized proportional to the mean expression.
\cr\cr
The proportion missing plot highlights samples that have elevated levels of missing.  Elevated missing could be indicative of a low RNA, inhibition or tissue type.  
\cr\cr
The rna content plot compares estimates of the total amount of RNA.  Any inconsistency between housekeeping and global RNA content estimates should be evaluated for assay issues.  In the case of miRNA's this could be due to the housekeeping genes undergoing different processing (no ligation).
\cr \cr
The batch effect plot looks for differences in sample summary statistics among groups, as defined by traits supplied during normalization.  Green dots are sized proportionally to their p-value and are less that 0.05.  A cartridge variable is automatically added to the checked traits.  If no header is supplied then the cartridge is assumed to be determined by column position in the input dataset.  That is, sample 1-12 from cartridge 1, 2-24 from cartridge 2 etc.  Note that the mean, sd and missing are post normalization and therefore may show an opposite trend than expected due to overcompensating normalization factors.  To check the raw data normalize using 'none' for each option.  Ideally, technical or study design traits should have orange dots and be near the line.  It can be expected that biological traits such as tumour status have different mean, sd and missing but they should not be correlated with assay related controls or normalization factors (positive, negative, housekeeping controls).  It's a good idea to double check any significant associations.
\cr \cr
The normalization factors plot looks at the calculated normalization parameters for each sample. By default samples are labeled if one of the normalization factors is greater or less than 100\% from the mean.  The raw data should be investigated for these samples as they could be contributing unnecessary noise to the downstream analysis.  Remember that the normalization factors are ratios of all samples vs a specifcic sample therefore high normalization factors reflect low code counts.  If a sample had a high positive control normalization factor this would show up as low intercept on the positive control plots.  The trend will also be opposite for the batch effects plots which use the summary of the counts not the ratio based normalization factor.  Normalization factors outside 100\% difference from the mean could reflect noise inducing over or under correction.
\cr \cr
The positive control plots compare the expected concentration with observed counts.  Orange dots are the negative controls.  The intercept reflects the sensitivity of the assay for that sample. 

}
\examples{

# load data
data(NanoString);

# specifiy housekeeping genes in annotation 
NanoString.mRNA[NanoString.mRNA$Name \%in\% 
	c('Eef1a1','Gapdh','Hprt1','Ppia','Sdha'),'Code.Class'] <- 'Housekeeping';

# setup the traits
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

# normalize
NanoString.mRNA.norm <- NanoStringNorm(
	x = NanoString.mRNA,
	anno = NA,
	CodeCount = 'geo.mean',
	Background = 'mean.2sd',
	SampleContent = 'housekeeping.geo.mean',
	round.values = TRUE,
	take.log = TRUE,
	traits = trait.strain,
	return.matrix.of.endogenous.probes = FALSE
	);

# plot all the plots as PDF report
pdf('NanoStringNorm_Example_Plots_All.pdf')
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = 'all'
	);
dev.off()

# publication quality tiff volcano plot
tiff('NanoStringNorm_Example_Plots_Volcano.tiff', units = 'in',  height = 6, 
	width = 6, compression = 'lzw', res = 1200, pointsize = 10);
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = c('volcano'),
	title = FALSE
	);
dev.off()

# all plots as seperate files output for a presentation
png('NanoStringNorm_Example_Plots_\%03d.png', units = 'in',  height = 6,
	width = 6, res = 250, pointsize = 10);
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = TRUE,
	plot.type = c('cv','mean.sd','RNA.estimates','volcano','missing','norm.factors','positive.controls','batch.effects')
	);
dev.off()

# user specified labelling with optimal resolution for most digital displays
png('NanoStringNorm_Example_Plots_Normalization_Factors.png', units = 'in', height = 6,
	width = 6, res = 250, pointsize = 10);
Plot.NanoStringNorm(
	x = NanoString.mRNA.norm,
	label.best.guess = FALSE,
	label.ids = list(genes = rownames(NanoString.mRNA.norm$gene.summary.stats.norm), 
		samples = rownames(NanoString.mRNA.norm$sample.summary.stats)),
	plot.type = c('norm.factors')
	);
dev.off()

}
\author{Daryl M. Waggott}

\keyword{NanoString}
\keyword{Normalization}
\keyword{mRNA}  
\keyword{miRNA}
