test.sample.summary.stats <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){
	
	# go to test data directory
	path.to.input.files <- '../NanoStringNorm/extdata/input_function_files/';
	path.to.output.files <- '../NanoStringNorm/extdata/output_function_files/';

	# read input files
	x             <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, '2011-10-01_NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.sample.summary <- dget(file = paste(path.to.output.files, date.checked.output, '_NanoString_mRNA_TCDD_Sample_Summary.txt', sep = ''));

	# run function to get *test output* 
	test.output.sample.summary      <- NanoStringNorm:::get.sample.summary.stats(x, anno);

	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.sample.summary, test.output.sample.summary);

	### check2 - check garbage input
	#check2.1 <- checkException(NanoStringNorm:::code.count.normalization(x, anno, 'mean'));

	checks <- c(check1.1 = check1.1);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
