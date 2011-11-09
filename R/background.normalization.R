background.normalization <- function(x, anno, Background = 'none', verbose = TRUE) {

	# check if missing
	if (is.na(Background)) {
		stop('Background: Background normalization method cannot be missing.  Try setting to *none*');		
		}

	# Background Correction
	if (Background != 'none') {

		# take the mean of negative controls
		if (Background == 'mean') {
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',], 
				MARGIN = 2, 
				FUN = mean, 
				na.rm = TRUE
				);
			}

		# take the mean of negative controls + 2 sd
		else if (Background == 'mean.2sd') {
			mean.plus.2sd <- function(y) mean(y, na.rm = TRUE) + 2 * sd(y, na.rm = TRUE);
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',], 
				MARGIN = 2, 
				FUN = mean.plus.2sd
				);
			}

		# take the maximum of negative controls
		else if (Background == 'max'){
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',], 
				MARGIN = 2, 
				FUN = max, 
				na.rm = TRUE
				);
			}

		# give an error if an unknown method is used
		else {
			stop('Background: Unimplemented Background method');
			}

		# flag any samples that are outliers
		background.level.sd.from.mean <- data.frame(background.zscore = (background.level - mean(background.level)) / sd(background.level));
		#colnames(background.level.sd.from.mean) <- 'background.zscore';

		if (verbose & any(abs(background.level.sd.from.mean) > 3)) {
			cat('Background: The following samples have an estimated bacground greater than 3 standard deviations from the mean.\n\n');
			print(signif(subset(background.level.sd.from.mean, abs(background.level.sd.from.mean) > 3),3));
			cat('\n');
			}

		# subtract backgound
		x <- t(apply(
			X = x, 
			MARGIN = 1, 
			FUN = '-', 
			background.level
			));

		# set negative values to zero
		x[x < 0] <- 0;

		# calculate count of missing values (DW missing could be moved to own function)
		genes.proportion.missing <- apply(
			X = x[grep('Endogenous', anno$Code.Class),],
			MARGIN = 1,
			FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y) }
			);

		samples.proportion.missing <- apply(
			X = x[grep('Endogenous', anno$Code.Class),],
			MARGIN = 2,
			FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y) }
			);

		samples.proportion.missing <- data.frame(row.names = colnames(x), proportion.missing = samples.proportion.missing);
		#colnames(samples.proportion.missing) <- 'proportion.missing';

		if (verbose == TRUE) {
			cat(paste('Background: After correction' , sum(samples.proportion.missing <= .90), 'samples and', sum(genes.proportion.missing <= .90), 'Endogenous genes have less than 90% missing. \n\n'));
			}

		if ( any(samples.proportion.missing > 0.9) ) {
			print(signif(subset(samples.proportion.missing, samples.proportion.missing > 0.9,),3));
			cat('\n');
			}

		}

	return(
		list(
			x = x,
			background.level = background.level
			)
		);
	}
