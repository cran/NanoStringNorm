read.markup.RCC <- function(rcc.path = ".", rcc.pattern = "*.RCC|*.rcc", exclude = NULL, include = NULL) {

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
		rcc.header <- read.table(paste(rcc.path, rcc.file, sep = "/"), nrows = 15, comment.char = "<", sep = ",", as.is = TRUE);
		rcc.data <- read.table(paste(rcc.path, rcc.file, sep = "/"), skip = 25, header = TRUE, comment.char = "<", sep = ",", as.is = TRUE);

		sample.name <- gsub(".RCC", "", gsub(" ", "_", rcc.file));
		colnames(rcc.header)[2] <- sample.name;
		colnames(rcc.data)[4] <- sample.name;

		if (count == 1) {
			rcc.header.merged <- rcc.header;
			rcc.data.merged <- rcc.data;
			}
		else {
			rcc.header.merged <- data.frame(rcc.header.merged, subset(rcc.header, select = 2));
			rcc.data.merged <- data.frame(rcc.data.merged, subset(rcc.data, select = 4));
			}

		count=count+1;
		}

	rcc.header.merged[3,1] <- "sample.id";
	rcc.header.merged[9,1] <- "lane.id";
	rownames(rcc.header.merged) <- rcc.header.merged[,1];
	rcc.header.merged <- rcc.header.merged[,-1];

	colnames(rcc.data.merged[1]) <- "Code.Count";

	x <- list(x = rcc.data.merged, header = rcc.header.merged);

	class(x) <- 'NanoString';
	return(x);
	}
