other.normalization.rank.normal <- function(x, anno, OtherNorm = 'none', verbose = TRUE, genes.to.fit = NA, genes.to.predict = NA) { 
	
	genes.to.fit <- if(is.na(genes.to.fit)) "all" else genes.to.fit;
	genes.to.predict <- if(is.na(genes.to.predict)) "all" else genes.to.predict
	
	# parse the genes.to.fit options
	if (any(genes.to.fit == 'all')) {
		genes.to.fit <- unique(anno$Code.Class);
		}
	if (any(genes.to.fit == 'controls')) {
		genes.to.fit <- c('Housekeeping','Negative','Positive','Control');
		}

	# set the genes not to be normalized to NA
	genes.to.fit = anno$Name %in% genes.to.fit | grepl(paste(genes.to.fit, collapse = "|"), anno$Code.Class, ignore.case = TRUE) ;
	
	#set values not to be fit to NA i.e. ignored
	x.fit <- x;
	x.fit[!genes.to.fit,] <- NA;

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		if (any(!genes.to.fit)) {
			cat("OtherNorm.rank.normal: The following genes will not be processed:\n");
			print(anno[!genes.to.fit,c("Code.Class","Name")]);
			}
		}
	
	# rank (forced normal) normalization
	x.fit[x.fit == 0] <- NA;

	x <- apply(
		X = x.fit,
		MARGIN = 2,
		FUN = function(y) qnorm (p = (rank(y) -.5)/length(y))
		);

	# add back the original counts to genes that should be ignored 
	x.fit <- data.frame(x.fit);
	if (any(!genes.to.fit)) {
		x.fit[!genes.to.fit, ] <- x[!genes.to.fit,];
		}

	rownames(x.fit) <- anno$Name;
	
	return(x.fit);
	}
