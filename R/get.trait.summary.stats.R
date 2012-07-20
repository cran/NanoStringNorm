get.trait.summary.stats <- function(x, anno = NA, logged, traits) {

	# attempt to convert traits into a matrix
	traits <- as.matrix(traits);

	# initialize values
	trait.summary.stats <- NULL;

	# loop over each trait and get the ttest pvalue and foldchange
	for ( i in 1:ncol(traits) ) {
		trait.value <- as.numeric(traits[,i]);
		trait.summary <- NULL;

		# group1 is reference (Normal) and group2 is effect (Tumour)
		trait.summary <- apply(
			X = x,
			MARGIN = 1,
			FUN = NanoStringNorm:::get.ttest.and.foldchange,
			group1 = (trait.value == 1),
			group2 = (trait.value == 2),
			logged = logged
			);

		# add the results to a matrix
		trait.summary <- t(trait.summary);
		trait.summary[is.nan(trait.summary)] <- NA;
		trait.summary.stats <- cbind(trait.summary.stats, trait.summary);

		}

	# round data to 2 sigfigs and add column nanes
	trait.summary.stats <- signif(trait.summary.stats, 2);
	colnames(trait.summary.stats) <- paste(c("P_", "FC_"), rep(colnames(traits), each = 2), sep = "");
	trait.summary.stats[is.nan(trait.summary.stats[,1]),c(1,2)] <- NA;

	return(trait.summary.stats);
	} 

