other.normalization <- function(x, anno, OtherNorm = 'none', verbose = TRUE, genes.to.fit = NA, genes.to.predict = NA, ...) { 

	if (OtherNorm == 'quantile') {
		x = NanoStringNorm:::other.normalization.quantile(x, anno, OtherNorm, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (OtherNorm == 'zscore') {
		x = NanoStringNorm:::other.normalization.zscore(x, anno, OtherNorm, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (OtherNorm == 'rank.normal') {
		x = NanoStringNorm:::other.normalization.rank.normal(x, anno, OtherNorm, verbose = verbose, genes.to.fit = genes.to.fit);
		}
	else if (OtherNorm == 'vsn') {
		x = NanoStringNorm:::other.normalization.vsn(x, anno, OtherNorm, verbose = verbose, genes.to.fit = genes.to.fit, genes.to.predict);
		}
	else if (OtherNorm == 'none') {
		return(x);
		}
	else {
		stop(paste('OtherNorm:  The OtherNorm option', OtherNorm, 'is not implemented try using one of quantile, zscore, rank.normal, vsn.'));
		}

	return(x);
	}
