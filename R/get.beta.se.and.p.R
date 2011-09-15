get.beta.se.and.p <- function(x, variable) {
	tryCatch(
		expr = c(	
			coef(summary(
				lm(variable ~ x)
				))[2,c(1,2,4)]
			),
		 error = function(e) { return ( c(NA, NA, NA)); }
		 );
	}


