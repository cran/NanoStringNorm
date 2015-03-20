 \name{Plot.NanoStringNorm.gvis}
 \title{Plot.NanoStringNorm.gvis}
 \alias{Plot.NanoStringNorm.gvis}
 \description{Plot interactive normalization diagnostic's and differential expression results using Google charts}
 \usage{
 Plot.NanoStringNorm.gvis(
         x,
         plot.type = c('gene.norm', 'sample'),
         save.plot = FALSE,
         path.to.mongoose = 'web',
         output.directory = 'NanoStringNorm_gvis_plots' 
         );
 
 }
 \arguments{
   \item{x}{An object of class NanoStringNorm i.e. the list output of NanoStringNorm}
   \item{plot.type}{A vector of the different plots that you would like.  Currently the following are implemented: gene.norm, gene.raw and sample.  These correspond to the list items in the NanoStringNorm     output. }
   \item{save.plot}{By default the plot is displayed on your browser using a built in web server.  This is good for quick evaluation if you have access to the underlying data.   However for distribution and presentations it is easier to save the plots bundled with an embedded web server.  You need to set save.plot to TRUE in order to get the bundled output.}
   \item{path.to.mongoose}{The path to a previously downloaded version of the mongoose embedded web server. \url{http://code.google.com/p/mongoose/} If the path to mongoose is not specified then it will be attempted to be downloaded from the internet i.e. default path.to.mongoose = 'web'. If you do not want to attempt to download from the internet then specify 'none'. }
   \item{output.directory}{The name of the directory you want the output saved in. }
}
\value{
 The output is a series of html files that contain interactive 4 dimensional plots.  The motivation for these plots is to enable end users to explore their own data quickly and effectively without depending on informaticians.
 \cr \cr
 For samples, a Google chart is rendered using sample summary statistics, traits and normalization factors as output from the NanoStringNorm function.  You can change the variables associated with each axis, colour coding and point sizing.  There are a number of ways to query the data in order to emphasise specific samples.  You can however mouse over the point to get label information or select specific rows from the table.
 \cr \cr
 Similarly, for gene level data you can explore the results summary statistics and trait associations.  For example, change the x-axis to Fold-Change.trait.X and the y-axis to P-value.trait.x, the size of the dot to mean expression and the colour to Code.Class.  This gives a nice interactive volcano plot.  Alternatively, you can plot the results of P-value.trait.x vs P-value-trait.y in order to look for common genes effecting different traits.  Note that P-values are displayed as -log10 P-value.
 \cr \cr
 To display plots saved as a bundle, first navigate to the directory and execute the mongoose binary.  Then start your browser and specify the following url \file{http://127.0.0.1:8080}  There should be a set of links to the plots.  The windows executable binary is automatically downloaded if no location is provided.  For other systems the source code is downloaded.  To install untar the binary, run make and then copy the resulting binary to the directory containing the output html files. 
 \cr \cr
 \code{cd output.directory;}
 \cr
 \code{tar -zxvf mongoose.tgz;}
 \cr
 \code{cd mongoose;}
 \cr
 \code{make linux;}
 \cr
 \code{cp mongoose ..;}
 \cr
 \code{./mongoose;}
 \cr \cr
 For further documentation regarding these plots see the R package googleVis and Google's visualization API.
 \cr \cr
 \url{http://code.google.com/p/google-motion-charts-with-r}
 \cr
 \url{http://code.google.com/apis/visualization/interactive_charts.html}
}
 
 \examples{
 
 # load data
 data(NanoString);
 
 # specify housekeeping genes in annotation 
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
 
 # plot the sample summaries to your browser
 if (requireNamespace("googleVis")) \{
 Plot.NanoStringNorm.gvis(
         x = NanoString.mRNA.norm,
         plot.type = c('gene.norm', 'sample'),
         save.plot = FALSE
         );
 
 # plot the gene summaries to a directory for distribution and later viewing
 Plot.NanoStringNorm.gvis(
        x = NanoString.mRNA.norm,
         plot.type = c('gene.norm', 'sample'),
         save.plot = TRUE,
         path.to.mongoose = 'none',
         output.directory = "NanoStringNorm_Interactive_Plot"
         );
 \}
 }
 \author{Daryl M. Waggott}
 
 \keyword{NanoString}
 \keyword{Normalization}
 \keyword{mRNA}
 \keyword{miRNA}



