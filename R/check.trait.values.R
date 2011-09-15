check.trait.values <- function(x, anno = NA, log, traits) {

	# attempt to convert traits into a matrix
	traits <- as.matrix(traits);

	# initialize variables
	check.na <- check.values <- check.samples <- 0;

	# run some checks
	if ( all(dim(traits) == 1) ) {
		# check single values specifically if NA
		check.na <- is.na(traits);
		}
	else {
		# check the data is numeric with only values of NA, 1 or 2
		check.values <- is.numeric(traits) & all(traits %in% c(NA,1,2));
		check.samples <- nrow(traits) == ncol(x) & all(rownames(traits) == colnames(x));
		}

	# evaluate if trait input was good
	if ( check.na ) {
		check <- 0;
		}
	else if ( check.values & check.samples ) {
		check <- 1;
		}
	else {
		# print warning messages
		if ( !check.na ) stop("Trait: Unrecognized trait input.");
		if ( !check.values ) stop("Trait: Only numeric variables with values NA, 1 or 2 are accepted.  The effect is terms of the second level i.e. disease.");
		if ( !check.samples ) stop("Trait: Values must include the same samples as the NanoString input data.  Confirm that traits have rownames in the same order as the columns in the expression data.");
		stop("Trait: Problem with trait data.  Check the documentation.");
		}

	return(check)
	}
