test.code.count.normalization <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){
	
	# go to test data directory
	path.to.input.files <- '../NanoStringNorm/extdata/input';
	path.to.output.files <- '../NanoStringNorm/extdata/output';

	# read input files
	x             <- read.table(paste(path.to.input.files, 'mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, 'mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, 'mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.predict.conc <- dget(file = paste(path.to.output.files, 'mRNA_TCDD_Predict_Conc.txt', sep = ''));
	
	# run function to get *test output* 
	test.output.predict.conc      <- NanoStringNorm:::predict.concentration(x, anno, log = TRUE, verbose = FALSE);

	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.predict.conc, test.output.predict.conc);

	### check2 - check garbage input
	#check2.1 <- checkException(NanoStringNorm::predict.concentratin(x, anno, 'mean', verbose = FALSE));

	checks <- c(check1.1 = check1.1);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
