code.count.normalization <- function(x, anno, CodeCount = 'none', verbose = TRUE) {

	# check if missing
	if (is.na(CodeCount)) {
		stop('CodeCount: CodeCount normalization method cannot be missing.  Try setting to *none*');		
		}

	# Code Count Normalization using sum: take sum of positive controls per sample
	if ( CodeCount == 'sum' ) { 
		pos.sample <- apply(
			X = x[anno$Code.Class == 'Positive',], 
			MARGIN = 2, 
			FUN = sum
			);
		}
	
	# Code Count Normalization using geometric mean
	else if ( CodeCount == 'geo.mean' ) { 
		pos.sample <- apply(
			X = x[anno$Code.Class == 'Positive',],
			MARGIN = 2, 
			FUN = NanoStringNorm:::get.geo.mean
			);
		}
		
	else {
		stop('CodeCount: Unimplemented CodeCount normalization method');
		}
	
	# calculate the normalization factor as the ratio of mean vs sample values
	pos.norm.factor <- mean(pos.sample) / pos.sample;

	# warning if norm factor is outside expected range
	if ( any(pos.norm.factor > 3) | any(pos.norm.factor < 0.3) ) {
	
		pos.norm.factor.flagged.samples <- as.data.frame(pos.norm.factor[(pos.norm.factor > 3) | (pos.norm.factor < 0.3)]);
		colnames(pos.norm.factor.flagged.samples) <- 'pos.norm.factor';

		if (verbose) {
			cat('CodeCount: The following samples have positive normalization factors outside the recommended range of (0.3 to 3).  Consider removing them.\n\n');
			print(signif(pos.norm.factor.flagged.samples,3));
			cat('\n');
			}
		}

	# multiply normalization factor to data values.  output needs to be transposed to get correct orientation
	x <- t(
		apply(
			X = x, 
			MARGIN = 1, 
			FUN = '*', 
			pos.norm.factor
			)
		);

	return(
		list(
			x = x,
			pos.norm.factor = pos.norm.factor,
			pos.sample = pos.sample
			)
		);
	}
