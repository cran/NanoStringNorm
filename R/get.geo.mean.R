get.geo.mean <- function(y) {
	y[y < 1] <- 1;
	exp(mean(log(na.omit(y))));
	}
