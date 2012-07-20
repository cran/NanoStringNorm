get.geo.mean <- function(y, logged = FALSE) {
	if (logged){
		mean(na.omit(y));
		}
	else {
		y[y < 1] <- 1;
		exp(mean(log(na.omit(y))));
		}
	}
