test.sample.content.normalization <- function(date.input = '2011-11-04', date.checked.output = '2011-11-04'){

	# data directories
	path.to.input.files <- '../NanoStringNorm/extdata/input_function_files/';
	path.to.output.files <- '../NanoStringNorm/extdata/output_function_files/';

	# read input files
	x             <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_matrix.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	anno          <- read.table(paste(path.to.input.files, date.input, '_NanoString_mRNA_TCDD_anno.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);
	trait         <- read.table(paste(path.to.input.files, '2011-10-01_NanoString_mRNA_TCDD_strain_info.txt', sep = ''), sep = '\t', header = TRUE, as.is = TRUE);

	# read *checked output*
	checked.output.hk.sum <- dget(file = paste(path.to.output.files, date.checked.output,'_NanoString_mRNA_TCDD_housekeeping.sum_SampleContent_Normalization.txt', sep = ''));
	checked.output.hk.geo.mean <- dget(file = paste(path.to.output.files, date.checked.output,'_NanoString_mRNA_TCDD_housekeeping.geo.mean_SampleContent_Normalization.txt', sep = ''));
	checked.output.top.geo.mean <- dget(file = paste(path.to.output.files, date.checked.output,'_NanoString_mRNA_TCDD_top.geo.mean_SampleContent_Normalization.txt', sep = ''));
	checked.output.top.mean <- dget(file = paste(path.to.output.files, date.checked.output,'_NanoString_mRNA_TCDD_top.mean_SampleContent_Normalization.txt', sep = ''));
	checked.output.total.sum <- dget(file = paste(path.to.output.files, date.checked.output,'_NanoString_mRNA_TCDD_total.sum_SampleContent_Normalization.txt', sep = ''));

	# run function to get *test output* 
	test.output.hk.sum      <- NanoStringNorm:::sample.content.normalization(x, anno, 'housekeeping.sum', verbose = FALSE);
	test.output.hk.geo.mean <- NanoStringNorm:::sample.content.normalization(x, anno, 'housekeeping.geo.mean', verbose = FALSE);
	test.output.top.geo.mean <- NanoStringNorm:::sample.content.normalization(x, anno, 'top.geo.mean', verbose = FALSE);
	test.output.top.mean <- NanoStringNorm:::sample.content.normalization(x, anno, 'top.mean', verbose = FALSE);
	test.output.total.sum <- NanoStringNorm:::sample.content.normalization(x, anno, 'total.sum', verbose = FALSE);
	
	### check1 - compare checked output == test output
	check1.1 <- checkEquals(checked.output.hk.sum, test.output.hk.sum);
	check1.2 <- checkEquals(checked.output.hk.geo.mean, test.output.hk.geo.mean);
	check1.3 <- checkEquals(checked.output.top.geo.mean, test.output.top.geo.mean);
	check1.4 <- checkEquals(checked.output.top.mean, test.output.top.mean);
	check1.5 <- checkEquals(checked.output.total.sum, test.output.total.sum);

	### check2 - test bad input 
	check2.1  <- checkException(NanoStringNorm:::sample.content.normalization(x, anno, 'garbage', verbose = FALSE));

	checks <- c(check1.1 = check1.1, check1.2 = check1.2, check1.3 = check1.3, check1.4 = check1.4,check1.5 = check1.5, check2.1 = check2.1);
	if (!all(checks)) print(checks[checks == FALSE]);

	return(all(checks))

	}
