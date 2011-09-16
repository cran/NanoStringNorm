check.trait.values <- function(x, anno = NA, log, traits) {

	# attempt to convert traits into a matrix
	traits <- as.matrix(traits);

	# initialize variables
	check.na <- check.values <- check.samples <- NA;

	# check if input is a single value specifically if NA
	if ( all(dim(traits) == 1) ) {
		if (is.na(traits)) {
			traits.ok <- FALSE;check.trait.values
			}
		else {
			stop("Trait: Unrecognized trait input.");
			}
		}
	else {
		# check if the traits contain only 1 and 2
		if ( is.numeric(traits) & all(traits %in% c(NA,1,2)) ) {
			traits.ok <- TRUE;
			}
		else {
			stop("Trait: Only numeric variables with values NA, 1 or 2 are accepted.  The effect is terms of the second level i.e. disease.");
			}
		# check if traits have the right ids in the right order
		if ( nrow(traits) == ncol(x) & all(rownames(traits) == colnames(x)) ) {
			traits.ok <- TRUE
			}
		else {
			stop("Trait: Values must include the same samples as the NanoString input data.  Confirm that traits have rownames in the same order as the columns in the expression data.")
			}

		# count the number of unique values per trait
		n.values <- apply(X = traits, MARGIN = 2, FUN = function(y) {length(unique(na.omit(y)))});

		if ( any(n.values == 1) ) {
			stop("Trait: Some of your traits only have one value.");
			}
		else {
			traits.ok <- TRUE;
			}

		}

	return(traits.ok)
	}
