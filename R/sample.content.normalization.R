sample.content.normalization <- function(x, anno, SampleContent = 'none', verbose = TRUE) {

	# check if missing
	if (is.na(SampleContent)) {
		stop('SampleContent: SampleContent normalization method cannot be missing.  Try setting to *none*');		
		}

	# Sample Content Normalization
	if (SampleContent != 'none') {

		# take the mean of the top 75-120 rnas.exclude viral genes from normalization factors due to potential associatio with phenotype
		endogenous.genes <- grepl('Endogenous', anno$Code.Class) & !grepl('bkv-|ebv-|hcmv-|hiv1-|hsv1-|hsv2-|kshv-|mcv-',anno$Name);

		# Check for Endogenous: take sum of housekeeping genes.
		if (SampleContent == 'housekeeping.sum') {
			rna.content <- apply(
				X = x[anno$Code.Class %in% c('Control', 'Housekeeping', 'housekeeping'), ], 
				MARGIN = 2, 
				FUN = sum
				);
			}
	
		else if (SampleContent == 'housekeeping.geo.mean') {

			# take the geometric mean of the housekeeping genes.
			rna.content <- apply(
				X = x[anno$Code.Class %in% c('Control', 'Housekeeping', 'housekeeping'), ], 
				MARGIN = 2, 
				FUN = NanoStringNorm:::get.geo.mean
				);
			}

		# take mean of all counts
		else if (SampleContent == 'total.sum') {
			rna.content <- apply(
				X = x[endogenous.genes,], 
				MARGIN = 2, 
				FUN = sum
				);
			}
	
		else if (SampleContent == 'top.mean') {

			# sum each RNA 
			sum.rna <- apply(
				X = x[endogenous.genes,],
				MARGIN = 1,
				FUN = sum
				);

			# get the ranks of the RNA sums
			rank.rna <- (length(sum.rna) + 1) - rank(sum.rna, ties = 'first');

			rna.content <- apply(
				X = x[endogenous.genes, ][rank.rna <= 75,],
				MARGIN = 2,
				FUN = mean,
				na.rm = TRUE
				);
			}

		else if (SampleContent == 'top.geo.mean') {

			# sum each RNA
			sum.rna <- apply(
				X = x[endogenous.genes,],
				MARGIN = 1,
				FUN = sum
				);

			# get the ranks of the RNA sums
			rank.rna <- (length(sum.rna) + 1) - rank(sum.rna, ties = 'first');

			rna.content <- apply(
				X = x[endogenous.genes,][rank.rna <= 75,],
				MARGIN = 2,
				FUN = NanoStringNorm:::get.geo.mean
				);
			}

		else {
			stop('SampleContent: Unimplemented SampleContent method');
			}

		# calc normalization factor
		sampleContent.norm.factor <- mean(rna.content) / rna.content;

		# flag any samples that are outliers
		rna.content.sd.from.mean <- data.frame(rna.zscore = (rna.content - mean(rna.content)) / sd(rna.content));

		if (verbose & any(abs(rna.content.sd.from.mean) > 3)) {
			cat('SampleContent: The following samples have a normalization factor greater than 3 standard deviations from the mean.\n\n');
			print(signif(subset(rna.content.sd.from.mean, abs(rna.content.sd.from.mean) > 3),3));
			cat('\n');
			}

		# adjust the data based on the normalization factor
		x <- t(apply(
			X = x,
			MARGIN = 1,
			FUN = '*',
			sampleContent.norm.factor
			));
		}
	return(
		list(
			x = x,
			sampleContent.norm.factor = sampleContent.norm.factor,
			rna.content = rna.content
			)
		);
	}
