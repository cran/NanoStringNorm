
get.mean.sd.and.cv <- function(x) { 
	mean.sd.and.cv <- c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100));
	mean.sd.and.cv[is.nan(mean.sd.and.cv)] <- NA;
	return( mean.sd.and.cv ); 
	}


