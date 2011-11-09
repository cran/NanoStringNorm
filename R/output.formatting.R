output.formatting <- function(x, anno, otherNorm = 'none', round.values = FALSE, log = FALSE, verbose = TRUE) {

	if ( (!is.logical(round.values) | !is.logical(log)) | is.na(round.values) | is.na(log) ) {
		stop('Formatting: Round and Log need to be TRUE or FALSE\n\n');
		}

	# round the data to discrete values.
	if (round.values == TRUE) {
		x <- round(x, digits = 0);
		}

	# shouldn't be log-transforming z-scores
	if (log == TRUE & otherNorm == 'zscore') {
		stop('Formatting: log-transformation of z-scores is intentionally not implemented.  Set log = FALSE.\n\n');
		}

	# log2 transformation
	if (log == TRUE & otherNorm != 'zscore') {

		# handle negative values
		if (any(x < 1)) {

			if (verbose) {
				cat('log: Setting values less than 1 to 1 in order to calculate the log in positive space.\n\n');
				}

			# set all values less than one on the raw scale to 0 on the log scale i.e. undetected
			x[x < 1] <- 1;
			}

		# log-transform (base2 by default)
		x <- signif(log2(x), 4);
		}

	return(x);
	}
