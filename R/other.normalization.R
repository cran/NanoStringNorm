other.normalization <- function(x, anno, otherNorm, verbose = TRUE) {

	# run additional normalization functions
	if (otherNorm != 'none') {

		# quantile normalization.  applied to all probes.
		if (otherNorm == 'quantile') {

			# zeros are considered missing so ignored for empirical distribution generation and application
			x.na <- x;
			x.na[x.na == 0] <- NA;

			# sort each sample independently into a new dataset
			x.sort <- apply(
				X = x.na, 
				MARGIN = 2,
				FUN = sort,
				na.last = TRUE
				);

			# take median across samples for each rank and use this as an empirical distribution.
			empirical.distribution <- apply(
				X = x.sort, 
				MARGIN = 1,
				FUN = median,
				na.rm = TRUE
				);

			# get the ranks in the unsorted data
			x.rank <- apply(
				X = x.na, 
				MARGIN = 2,
				FUN = rank,
				na.last = TRUE,
				ties = 'first'
				);

			# set ranks to NA if NA in the original data
			x.rank[is.na(x.na)] <- NA;

			# convert the ranks into probabilities to use quantile functio
			get.prob.from.rank <- function(xval) ( xval - 1 ) / (length(na.omit(xval)) - 1);

			x.rank.prob <- apply(
				X = x.rank,
				MARGIN = 2,
				get.prob.from.rank
				);

			# use the empirical distribution to peg the ranks too
			get.quantile <- function(y, dist) { quantile(x = dist, probs = y, na.rm = TRUE) }

			# only run the quantile function on the endogenous but add back the controls
			x <- apply(
				X = x.rank.prob,
				MARGIN = 2, 
				FUN = get.quantile, 
				dist = empirical.distribution
				);

			# NA's back to 0
			x[is.na(x)] <- 0;
			}

		# z-score normalization
		else if (otherNorm == 'zscore') {
			x <- apply(
				X = x,
				MARGIN = 2,
				FUN = scale
				);
			}
		# give an error if the method is not known
		else {
			stop('Unimplemented otherNorm method');
			}
		}

	return(x);
	}
