get.gene.summary.stats <- function(x, anno = NA) {

	gene.summary.stats <- apply(
		X = x,
		MARGIN = 1,
		FUN = NanoStringNorm:::get.mean.sd.and.cv
		);

	gene.summary.stats <- t(gene.summary.stats);

	# calculate count of missing values (DW missing could be moved to own function)
	genes.proportion.missing <- apply(
		X = x,
		MARGIN = 1,
		FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y) }
		);

	gene.summary.stats <- cbind(gene.summary.stats, 100 * genes.proportion.missing);

	gene.summary.stats <- round(gene.summary.stats, 1);
	colnames(gene.summary.stats) <- c("Mean", "SD", "CV", "Missing");
	rownames(gene.summary.stats) <- anno$Name;

	return(gene.summary.stats);
	}


