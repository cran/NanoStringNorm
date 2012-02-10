other.normalization <- function(x, anno, otherNorm = 'none', verbose = TRUE, genes.to.fit, genes.to.predict, ...) { 
	
	if (otherNorm == 'quantile') {
		x = other.normalization.quantile(x, anno, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (otherNorm == 'zscore') {
		x = other.normalization.zscore(x, anno, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (otherNorm == 'rank.normal') {
		x = other.normalization.rank.normal(x, anno, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (otherNorm == 'vsn') {
		x = other.normalization.vsn(x, anno, verbose = verbose, genes.to.fit = genes.to.fit, genes.to.predict = genes.to.predict);
		}
	else {
		stop(paste('otherNorm:  The otherNorm option', otherNorm, 'is not implemented try using one of quantile, zscore, rank.normal, vsn.'));
		}

	#rownames(x) <- anno$Name;
	return(x);
	}
