get.trait.summary.stats <- function(x, anno = NA, logged, traits) {

	# attempt to convert traits into a matrix
	traits <- as.data.frame(traits,as.is=TRUE);

	# initialize values
	trait.summary.stats <- NULL;
	if ('pair.ids' %in% colnames(traits) && all(table(traits[,'pair.ids'])== 2)) {
#		pair.ids <- traits[,'pair.ids'];
		traits <- traits[order( traits[,'pair.ids']),];
		x <- x[,rownames(traits)];
		paired <- TRUE;
		traits <- traits[,!colnames(traits) %in% 'pair.ids', drop=FALSE];
		}
	else {
		paired <- FALSE;
		}

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
			logged = logged,
			paired = paired
			);

		# add the results to a matrix
		trait.summary <- t(trait.summary);
		trait.summary[is.nan(trait.summary)] <- NA;

		# don't forget to get the qvalues
		trait.summary <- cbind(trait.summary, p.adjust(trait.summary[,1], method = 'fdr'));

		# add the the results to a matrix
		trait.summary.stats <- cbind(trait.summary.stats, trait.summary);

		}

	# round data to 2 sigfigs and add column nanes
	trait.summary.stats <- signif(trait.summary.stats, 2);
	colnames(trait.summary.stats) <- paste(c("P_", "FC_", "Q_"), rep(colnames(traits), each = 3), sep = "");
	trait.summary.stats[is.nan(trait.summary.stats[,1]),c(1,2)] <- NA;

	return(trait.summary.stats);
	} 

