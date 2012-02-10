Plot.NanoStringNorm.gvis <- function(x, plot.type = c("gene.norm","gene.raw", "sample"), save.plot = FALSE, path.to.mongoose = NA, output.directory = "NanoStringNorm_gvis_plots") {


	if (!suppressPackageStartupMessages(require(googleVis))) {
		stop ("Plot.NanoStringNorm.gvis:  googleVis is not available.");
		}

	for (plot.item in plot.type) {

		if (plot.item == "gene.norm") {
			NSN.output.name <- "gene.summary.stats.norm";
			idvar <- "Gene";
			}
		else if (plot.item == "gene.raw") {
			NSN.output.name <- "gene.summary.stats.raw";
			idvar <- "Gene";
			}
		else if (plot.item == "sample") {
			NSN.output.name <- "sample.summary.stats.norm";
			idvar <- "Sample";
			}
		else {
			stop(paste("Plot.NanoStringNorm.gvis:  Unrecognized plot.type", plot.item))
			};

		data.to.plot <- x[[NSN.output.name]];
		
		# add the annotation and dummy time variable for classification in plotting
		if (grepl("gene", plot.item) ) { 
			data.to.plot <- cbind(x$raw.data[,c('Name','Code.Class')],data.to.plot);
			colnames(data.to.plot)[1] <- 'Gene';
			}
		else if (grepl("sample", plot.item)) {
			data.to.plot <- merge(data.to.plot, x$traits, by.x = 0, by.y = 0);
			colnames(data.to.plot)[1] <-  'Sample';
			# add a prefix to the sample names because they sometimes cause errors
			data.to.plot$Sample <- paste(1:nrow(data.to.plot), data.to.plot$Sample, sep = "-" );
			}

		data.to.plot$time <- 1;
		
		# take the -log10P for scaling.  This should only be done on gene.norm
		pval.columns <- colnames(data.to.plot)[grepl("P_",colnames(data.to.plot))];
		for (pval in pval.columns) {
			data.to.plot[,pval] <- round(-log10(data.to.plot[,pval]),2);
			}

		# the intial plotting parameters
		initial.plotting.parameters <- '{"iconKeySettings":[],"stateVersion":3,"time":"notime","xAxisOption":"_NOTHING","playDuration":15,"iconType":"BUBBLE","sizeOption":"_NOTHING","xZoomedDataMin":null,"xZoomedIn":false,"duration":{"multiplier":1,"timeUnit":"none"},"yZoomedDataMin":null,"xLambda":1,"colorOption":"_NOTHING","nonSelectedAlpha":0.4,"dimensions":{"iconDimensions":[]},"yZoomedIn":false,"yAxisOption":"_NOTHING","yLambda":1,"yZoomedDataMax":null,"showTrails":false,"xZoomedDataMax":null};';
		# "iconKeySettings":[{"key":{"dim0":"Cyp1b1"},"trailStart":"1901"}]

		# call googlevis and make motionChart
		plot.motion = gvisMotionChart(
			data = data.to.plot, 
			idvar = idvar, 
			timevar = "time", 
			options = list(
				gvis.editor="Editor",
				height=700,
				width=900,
				showChartButtons = TRUE,
				showHeader = TRUE,
				showSelectComponent = TRUE,
				showSidePanel = TRUE,
				showMetrixPicker = TRUE,
				showYMetricPicker = TRUE,
				showXScalePicker = TRUE,
				showYScalePicker = TRUE,
				showAdvancedPanel = TRUE,
				state = initial.plotting.parameters
				)
			);
		
		# create a data table containing the plotted data
		plot.table <- gvisTable(
			data = data.to.plot,
			options=list(width=600, height=700)
			);

		# merge the chart and table
		plot.merge <- gvisMerge(plot.motion, plot.table, horizontal = TRUE);

		# plot using internal R webserver
		plot(plot.merge);
		
		if (save.plot == TRUE) {
			
			# create a directory to dump the html files
			if (!file.exists(output.directory))  {
				dir.create(output.directory);
				}	
			
			# copy the mongoose embedded web server.  
			# note: mongoose was written under MIT licence http://code.google.com/p/mongoose/
			# if no path to mongoose then download the mongoose executables and put them in the report directory
			
			if (is.na(path.to.mongoose)) {
				download.file(url = "http://mongoose.googlecode.com/files/mongoose-3.0.exe", destfile = paste(output.directory,"/mongoose.exe",sep = ""));
				download.file(url = "http://mongoose.googlecode.com/files/mongoose-3.0.tgz", destfile = paste(output.directory, "/mongoose.tgz", sep = ""));
				cat("Plot.NanoStringNorm.gvis: Note that only the source code was downloaded for non windows systems.  You will have to untar and compile the code.  See the docs.\n");
				}
			else if (path.to.mongoose != "none") {
				if (file.exists(paste(path.to.mongoose, "/mongoose", sep = ""))) { 
					file.copy(from = paste(path.to.mongoose, "/mongoose"), to = output.directory);
					}
				if (file.exists(paste(path.to.mongoose, "/mongoose.exe", sep = ""))) { 
					file.copy(from = paste(path.to.mongoose, "/mongoose.exe"), to = output.directory);
					}
				}

			# Create Google Gadget
			#cat(createGoogleGadget(plot.merge), file = paste("NanoStringNorm_gvis_", plot.type ,"_summary.html", sep = ""))
			
			cat(plot.merge$html$chart, file=paste(output.directory, "/NanoStringNorm_gvis_", plot.item ,"_summary.html", sep = ""));

			# print message about links
			cat("Plot.NanoStringNorm.gvis: First run the mongoose binary found in the NanoStringNorm_gvis_plots and then navigate to http://127.0.01:8080 in your browser to view the plots\n");
			}

		}
	}
