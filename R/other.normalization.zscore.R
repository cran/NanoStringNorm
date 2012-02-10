other.normalization.zscore <- function(x, anno, otherNorm = 'none', verbose = TRUE, genes.to.fit = NULL, genes.to.predict = NULL) { 
	
	genes.to.fit <- if(is.null(genes.to.fit)) "endogenous" else genes.to.fit;
	genes.to.predict = if(is.null(genes.to.predict)) "endogenous" else genes.to.predict;

	# parse the genes.to.fit options
	if (any(genes.to.fit == 'all')) {
		genes.to.fit <- unique(anno$Code.Class);
		}
	if (any(genes.to.fit == 'controls')) {
		genes.to.fit <- c('Housekeeping','Negative','Positive','Control');
		}

	# which genes need to be fit
	genes.to.fit = anno$Name %in% genes.to.fit | grepl(paste(genes.to.fit, collapse = "|"), anno$Code.Class, ignore.case = TRUE) ;
	
	# set values not to be fit to NA i.e. ignored
	x.fit <- x;
	x.fit[!genes.to.fit,] <- NA;

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		if (any(!genes.to.fit)) {
			cat("OtherNorm.zscore: The following genes will not be processed:\n");
			print(anno[!genes.to.fit,c("Code.Class","Name")]);
			}
		}
	
	#### START METHOD #########################################################################
	# z-score (standard normal) normalization
	x.fit <- apply(
		X = x.fit,
		MARGIN = 2,
		FUN = scale
		);
	
	### END METHOD ############################################################################

	# add back the original counts to genes that should be ignored
	x.fit <- data.frame(x.fit);
	if (any(!genes.to.fit)) {
		x.fit[!genes.to.fit, ] <- x[!genes.to.fit,];
		}

	rownames(x.fit) <- anno$Name;

	return(x.fit);
	}
