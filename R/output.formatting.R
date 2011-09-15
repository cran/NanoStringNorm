output.formatting <- function(x, anno, otherNorm, round.values, log, verbose = TRUE) {
	
	# round the data to discrete values.
	if (round.values == TRUE) {
		x <- round(x, digits = 0);
		}

	# shouldn't be log-transforming z-scores
	if (verbose & log == TRUE & otherNorm == 'zscore') {
		cat('log: log-transformation of z-scores is intentionally not implemented\n\n');
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
		x <- log2(x);
		}

	return(x);
	}
