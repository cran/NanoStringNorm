get.sample.summary.stats <- function(x, anno) {

		# calculate mean of each sampleContent
		sample.mean <- apply(
			X = x[grep('Endogenous', anno$Code.Class),],
			MARGIN = 2,
			FUN = mean,
			na.rm = TRUE
			);

		# calculate SD of each sampleContent
		sample.sd <- apply(
			X = x[grep('Endogenous', anno$Code.Class),],
			MARGIN = 2,
			FUN = sd,
			na.rm = TRUE
			);

		# calculate proportion missing values
		sample.proportion.missing <- apply(
			X = x[grep('Endogenous', anno$Code.Class),],
			MARGIN = 2,
			FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y); }
			);

	return(
		data.frame(
			row.names = names(x),
			Sample.Mean = sample.mean,
			Sample.SD = sample.sd,
			Sample.Missing = sample.proportion.missing
			)
		);
	}
