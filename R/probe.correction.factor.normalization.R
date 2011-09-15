probe.correction.factor <- function(x, anno, Probe.Correction.Factor, verbose = TRUE) {

	# first convert it to a matrix
	Probe.Correction.Factor <- as.matrix(Probe.Correction.Factor);

	if ( !(all(Probe.Correction.Factor %in% c('none', 'filter')) | dim(Probe.Correction.Factor)[2] == 2) ) {
		stop("Probe.Correction: The parameter argument is not expected, see the documentation.");
		}

	# check for flags indicating required probe level background correction
	if ( any(grepl('Message', anno$Name, fixed = TRUE)) & all(Probe.Correction.Factor == 'none') ) {
		if (verbose) {
			cat('Genes requiring probe level correction:\n\n');
			print(as.character(anno[grepl('Message', anno$Name), 'Name']));
			cat('\n');
			}
		stop('Probe.Corection: Your data needs probe level background correction.  See documentation regarding Probe.Correction.Factor parameter.');
		}

	# filter genes with probe level warnings.  this can be use if there is no file for probe correction.
	if ( all(Probe.Correction.Factor == 'filter') ) {
		if (verbose) {
			cat('You are removing the following probes from all further analysis:\n\n');
			print(as.character(anno[grepl('Message', anno$Name), 'Name']));
			cat('\n');
			}
		x <- x[!grepl('Message', anno$Name, fixed = TRUE),];
		anno <- anno[!grepl('Message', anno$Name, fixed = TRUE),];
		}

	# return results if no probe correction
	if ( all(Probe.Correction.Factor == 'none') | all(Probe.Correction.Factor == 'filter') | length(Probe.Correction.Factor) == 1 ) {
		return(
			list(
				x = x,
				anno = anno
				)
			);
		}

	# first convert it to a data.frame
	Probe.Correction.Factor <- as.data.frame(Probe.Correction.Factor, stringsAsFactors = FALSE);

	# check validity of correction data
	if ( ncol(Probe.Correction.Factor) != 2 ) {
		stop('Probe.Correction:  There must be two columns including gene name and correction factor.  See documentation.');
		}
	
	# correct for probe level background
	if ( ncol(Probe.Correction.Factor) == 2 ) {
		colnames(Probe.Correction.Factor) <- c('Name', 'Probe.Correction.Factor');

		# remove probe level background flag from the gene names in order to check for matching ids
		Probe.Correction.Factor$Name <- gsub(' ', '', Probe.Correction.Factor$Name, fixed = TRUE);
		anno$Name <- gsub(' \\(.+$', '', anno$Name);
		anno$Name <- gsub(' ', '', anno$Name, fixed = TRUE);
		rownames(x) <- anno$Name;
		#anno <- anno[anno$Code.Class != 'Message',];
		#x <- x[anno$Code.Class != 'Message',];

		# check that the gene names match
		if ( !all(Probe.Correction.Factor$Name %in% anno$Name) ) {
			print(Probe.Correction.Factor[rownames(Probe.Correction.Factor) %in% anno$Name,]);
			stop('Probe.Correction: Not all probes in Probe.Correction.Factor are found in your data. Check the names.');
			}

		# expand background correction factor to all genes.  i.e. set other genes to 0
		Probe.Correction.Factor <- merge(anno[,'Name'], Probe.Correction.Factor, by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE, sort = FALSE);
		Probe.Correction.Factor[is.na(Probe.Correction.Factor$Probe.Correction.Factor),'Probe.Correction.Factor'] <- 0;

		# sort according to original order
		rownames(Probe.Correction.Factor) <- Probe.Correction.Factor[,1];
		Probe.Correction.Factor <- Probe.Correction.Factor[anno$Name,];
		Probe.Correction.Factor[,1] <- NULL;
		Probe.Correction.Factor <- as.numeric(Probe.Correction.Factor[,1]);

		# adjust for probe level background correction
		Probe.Correction.Factor <- sapply(X = x[anno$Name == 'POS_A(128)',], FUN = '*', Probe.Correction.Factor);
		x <- x - as.matrix(Probe.Correction.Factor);
		x[x < 0] <- 0;
		}
		
	return(
		list(
			x = x,
			anno = anno
			)
		);

	}
