# The NanoStringNorm package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

read.markup.RCC <- function(rcc.path = ".", rcc.pattern = "*.RCC|*.rcc", exclude = NULL, include = NULL, nprobes = -1) {

	# check if path exist
	if (! file.exists(rcc.path) ) {
		stop(paste("READ.MARKUP.RCC: path was not found.  \n")) ;
		}

	rcc.files <- list.files(path = rcc.path, pattern = rcc.pattern);
	rcc.files <- rcc.files[!rcc.files %in% exclude | rcc.files %in% include]

	# check if any rcc.files
	if (length(rcc.files) == 0) {
		stop(paste("READ.MARKUP.RCC: no RCC files found.  \n")) ;
		}

	rcc.header.merged <- NULL;
	rcc.data.merged <- NULL;

	count = 1;
	for (rcc.file in rcc.files) {

		cat("\nreading RCC file [", count, "/", length(rcc.files), "]: ", rcc.file, sep = "");

		# read RCC file and enclose in valid document tags for XML parser
		data <- xmlParse(
			paste(
				"<doc>", 
				paste(readLines(paste(rcc.path, rcc.file, sep = "/"), warn = FALSE), collapse = "\r"), 
				"</doc>", 
				sep = ""
				)
			);

		# extract info for various tags
		header.info <- xmlToDataFrame(getNodeSet(data, "/doc//Header"));
		sample.info <- xmlToDataFrame(getNodeSet(data, "/doc//Sample_Attributes"));
		lane.info <- xmlToDataFrame(getNodeSet(data, "/doc//Lane_Attributes"));
		code.info <- xmlToDataFrame(getNodeSet(data, "/doc//Code_Summary"));

		# convert data into table for all data structures
		header.con <- textConnection(as.character(header.info$text));
		header.data <- read.csv(header.con, header = FALSE, stringsAsFactors = FALSE);

		sample.con <- textConnection(as.character(sample.info$text));
		sample.data <- read.csv(sample.con, header = FALSE, stringsAsFactors = FALSE);
		sample.data[which(sample.data[, 1] == "ID"), 1] <- "sample.id";

		lane.con <- textConnection(as.character(lane.info$text));
		lane.data <- read.csv(lane.con, header = FALSE, stringsAsFactors = FALSE);
		lane.data[which(lane.data[, 1] == "ID"), 1] <- "lane.id";

		code.con <- textConnection(as.character(code.info$text));
		code.data <- read.csv(code.con, header = TRUE, stringsAsFactors = FALSE);

		# combine metadata
		rcc.header <- do.call(rbind, list(header.data, sample.data, lane.data));
		rcc.data <- code.data;

		# assign rownames for lookup
		rownames(rcc.header) <- rcc.header[, 1];

		# assign sample names
		sample.name <- gsub(".RCC", "", gsub(" ", "_", rcc.file));
		colnames(rcc.header)[2] <- sample.name;
		colnames(rcc.data)[4] <- sample.name;

		if (count == 1) {
			rcc.header.merged <- rcc.header;
			rcc.data.merged <- rcc.data;
			}
		else {
			rcc.header.merged <- data.frame(
				rcc.header.merged, 
				subset(rcc.header[rcc.header.merged[, 1], ], select = 2)
				);
			rcc.data.merged <- data.frame(rcc.data.merged, subset(rcc.data, select = 4));
			}

		count=count+1;
		}

	# assign rownames
	rownames(rcc.header.merged) <- rcc.header.merged[,1];
	rcc.header.merged <- rcc.header.merged[,-1];

	# set column name
	colnames(rcc.data.merged[1]) <- "Code.Count";

	# merge data structures
	x <- list(x = rcc.data.merged, header = rcc.header.merged);

	class(x) <- 'NanoString';

	return(x);

	}
