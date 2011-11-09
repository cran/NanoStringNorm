test.output.formatting <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){
	
	# go to test data directory
	path.to.input.files <- '../NanoStringNorm/extdata/input_function_files/';
	path.to.output.files <- '../NanoStringNorm/extdata/output_function_files/';

	# read input files
	x     <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno  <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait <- read.table(paste(path.to.input.files, '2011-10-01_NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.roundT.logT <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_Output_Formatting_roundT_logT.txt', sep = ''));

	# run function to get *test output* 
	test.output.roundT.logT <- NanoStringNorm:::output.formatting(x, anno, round.values = TRUE, log = TRUE, verbose = FALSE);

	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.roundT.logT, test.output.roundT.logT);

	### check2 - check garbage input
	check2.1 <- checkException(NanoStringNorm:::output.formatting(x, anno, round.values = 'garbage', verbose = FALSE));
	check2.2 <- checkException(NanoStringNorm:::output.formatting(x, anno, round.values = NA, verbose = FALSE));
	
	### check3 - check if log complains with zscore normalization (i.e. negatives)
	check3.1 <- checkException(NanoStringNorm:::output.formatting(x, anno, round.values = TRUE, log = TRUE, otherNorm = 'zscore', verbose = FALSE));

	checks <- c(check1.1 = check1.1, check2.1 = check2.1, check2.2 = check2.2);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
