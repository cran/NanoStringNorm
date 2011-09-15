predict.concentration <- function(x, anno, log, verbose = TRUE) {

	# concentration in fM
	if ( log == TRUE ) { 
		positive.controls.concentration <- c(.125, .5, 2, 8, 32, 128);
		}
	else {
		positive.controls.concentration <- log2(c(.125, .5, 2, 8, 32, 128));
		}

	# generate average model based on mean of postive controls
	endogenous.gene.mean <- apply(
		X = x[grepl('Endogenous', anno$Code.Class),],
		MARGIN = 1,
		FUN = mean,
		na.rm = TRUE
		);

	positive.control.mean <- apply(
		X = x[anno$Code.Class == 'Positive', ],
		MARGIN = 1,
		FUN = mean,
		na.rm = TRUE
		);

	fm1.intercept <- lm(positive.control.mean ~ positive.controls.concentration)$coefficients[1];
	fm1.slope     <- lm(positive.control.mean ~ positive.controls.concentration)$coefficients[2];

	mean.endogenous.gene.concentration <- (endogenous.gene.mean - fm1.intercept) / fm1.slope;

	# generate individual specific model of postive controls
	all.endogenous.gene.concentration <- NULL;
	for ( i in 1:ncol(x) ) {
		sample <- x[,i];
		fm2.intercept                     <- lm(sample[anno$Code.Class == 'Positive'] ~ positive.controls.concentration)$coefficients[1];
		fm2.slope                         <- lm(sample[anno$Code.Class == 'Positive'] ~ positive.controls.concentration)$coefficients[2];
		endogenous.gene.concentration     <- (sample[grepl('Endogenous', anno$Code.Class)] - fm2.intercept) / fm2.slope;
		all.endogenous.gene.concentration <- cbind(all.endogenous.gene.concentration, endogenous.gene.concentration);
		}

	colnames(all.endogenous.gene.concentration) <- colnames(x);

	gene.concentration <- cbind(
		row.names = anno[grepl('Endogenous', anno$Code.Class),'Name'], 
		mean.concentration = mean.endogenous.gene.concentration, 
		all.endogenous.gene.concentration
		)

	return(gene.concentration);
	}
