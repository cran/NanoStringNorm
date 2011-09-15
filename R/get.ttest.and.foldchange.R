get.ttest.and.foldchange <- function(x, group1, group2, log = TRUE) {
	group1[is.na(group1)] <- FALSE;
	group2[is.na(group2)] <- FALSE;
	
	tryCatch(
		expr = c(
			t.test(
				x[group1],
				x[group2],
				paired = FALSE,
				var.equal = FALSE,
				alternative = 'two.sided'
				)$p.value,
			if (log == TRUE) { round(mean(x[group2], na.rm = TRUE) - mean(x[group1], na.rm = TRUE),1); } else { round(mean(x[group2], na.rm = TRUE) / mean(x[group1], na.rm = TRUE),1); }
			),
		 error = function(e) { return ( c(NA, NA)); }
		 );
	}


