NanoStringNorm <- function(x, anno = NA, Probe.Correction.Factor = 'none', CodeCount = 'none', Background = 'none', SampleContent = 'none', otherNorm = 'none', round.values = FALSE, log = FALSE, return.matrix.of.endogenous.probes = FALSE, traits = NA, verbose = TRUE) {

	# format the annotation
	if ( all(is.na(anno)) ) {

		# check that the annotation columns exist
		if ( colnames(x)[1] == "CodeClass" ) colnames(x)[1] <- "Code.Class";

		if ( any(!c('Code.Class','Name', 'Accession') %in% colnames(x)) ) {
			stop ("You have not specified an annotation file and your data does not contain the Code.Class, Accession or Name fields.");
			}

		# remove probe level warning message line from data
		if ( grepl("+++ Functional tests", x[nrow(x),"Name"], fixed = TRUE) ){
			x <- x[-nrow(x),];
			}

		x.raw <- data.frame(x);

		anno <- x[,c('Code.Class', 'Name', 'Accession')];
		x <- as.matrix(x[,!names(x) %in% c('Code.Class', 'Name', 'Accession')]);
		rownames(x) <- anno$Name;

		}

	# start printing analysis log
	if (verbose) {
		cat('\n##############################\n');
		cat(paste('### NanoStringNorm v', packageDescription("NanoStringNorm")$Version, ' ###\n', sep = ''));
		cat('##############################\n\n');
		cat(paste('There are', ncol(x),'samples and', sum(grepl('Endogenous', anno$Code.Class)), 'Endogenous genes \n\n')); 
		}

	# Check data is a numeric matrix
	if (any(!is.numeric(x))) {
		stop('There are character values in your data.  Check that you labelled the annotation columns correctly');
		}

	# Check data in positive
	if (any(x < 0)) {
		stop('There are negative values in your data');
		}

	# Check for Endogenous probes
	if ( length(grep('Endogenous', anno$Code.Class)) == 0 ) {
		stop('There are no Endogenous genes in your data');
		}

	# check for Positive Controls
	if ( CodeCount != 'none' & !any(anno$Code.Class == 'Positive') ) {
		stop('You cannot do CodeCount Normalization. There are no Positive Controls in your data.');
		}

	# check Housekeeping and Control Genes
	if ( grepl('housekeeping', SampleContent) & !any(anno$Code.Class %in% c('Control', 'Housekeeping', 'housekeeping')) ) {
		stop('You Cannot do SampleContent Normalization. There are no *annotated* Housekeeping / Control genes in your data.');
		}

	# Check normalization parameters
	if ( !CodeCount %in% c('none', 'sum', 'geo.mean') ) {
		stop('Unrecognized CodeCount Normalization method');
		}
	if ( !Background %in% c('none', 'mean', 'mean.2sd', 'max') ) {
		stop('Unrecognized Background Normalization method');
		}
	if ( !SampleContent %in% c('none', 'housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'top.mean', 'top.geo.mean') ) {
		stop('Unrecognized SampleContent Normalization method');
		}
	if ( !otherNorm %in% c('none', 'quantile', 'zscore') ) {
		stop('Unrecognized otherNorm Normalization method');
		}

	# do Probe Correction Factor Normalization
	output.probe.correction.factor <- NanoStringNorm:::probe.correction.factor(x, anno, Probe.Correction.Factor, verbose);
	anno <- output.probe.correction.factor$anno;
	x <- output.probe.correction.factor$x;
	rm(output.probe.correction.factor);

	# get gene summary stats for raw data
	gene.summary.stats.raw  <- NanoStringNorm:::get.gene.summary.stats(x, anno);

	# do CodeCount Normalization
	if ( CodeCount %in% c('sum', 'geo.mean') ) {
		output.code.count.normalization <- NanoStringNorm:::code.count.normalization(x, anno, CodeCount, verbose);
		x <- output.code.count.normalization$x;
		pos.norm.factor <- output.code.count.normalization$pos.norm.factor;
		pos.sample <- output.code.count.normalization$pos.sample;
		rm(output.code.count.normalization);
		}

	# do Background Correction Normalization
	if ( Background %in% c('mean', 'mean.2sd', 'max') ) {
		output.background.normalization <- NanoStringNorm:::background.normalization(x, anno, Background, verbose);
		x <- output.background.normalization$x;
		background.level <- output.background.normalization$background.level;
		rm(output.background.normalization);
		}

	# do Sample Content Normalization
	if ( SampleContent %in% c('housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'top.mean', 'top.geo.mean') ) {
		output.sample.content.normalization <- NanoStringNorm:::sample.content.normalization(x, anno, SampleContent, verbose);
		x <- output.sample.content.normalization$x;
		sampleContent.norm.factor <- output.sample.content.normalization$sampleContent.norm.factor;
		rna.content <- output.sample.content.normalization$rna.content;
		rm(output.sample.content.normalization);
		}

	# do other additional normalizations.  note these are applied to all probes but excluding counts equal to 0
	if ( otherNorm %in% c('quantile', 'zscore') ) {
		x <- rbind(
			x[!grepl('Endogenous', anno$Code.Class),],
			NanoStringNorm:::other.normalization(x[grepl('Endogenous', anno$Code.Class),], anno, otherNorm, verbose)
			);
		}
	
	# do rounding, log-transformation
	x <- NanoStringNorm:::output.formatting(x, anno, otherNorm, round.values, log, verbose);

	# get predicted concentration based on positive controls
	predicted.concentration <- NanoStringNorm:::predict.concentration(x, anno, log, verbose);

	# output the data as a matrix filtering the annotation and control genes.
	if ( return.matrix.of.endogenous.probes == TRUE ) {
		return(x[grepl('Endogenous', anno$Code.Class),]);
		}

	# get sample summary stats
	sample.summary.stats <- NanoStringNorm:::get.sample.summary.stats(x, anno); 

	# add the normalization factors to the sample summary stats
	sample.summary.stats <- cbind(
		sample.summary.stats,
		pos.norm.factor = if (exists('pos.norm.factor')) pos.norm.factor else NA,
		pos.controls = if (exists('pos.sample')) pos.sample else NA,
		background.level = if (exists('background.level')) background.level else NA,
		sampleContent.norm.factor = if (exists('sampleContent.norm.factor')) sampleContent.norm.factor else NA,
		rna.content = if (exists('rna.content')) rna.content else NA,
		row.names = colnames(x)
		);

	# get gene summary stats for normalized data
	gene.summary.stats.norm <- NanoStringNorm:::get.gene.summary.stats(x, anno);

	# check that trait data is the right format
	check.traits <- NanoStringNorm:::check.trait.values(x, anno, log, traits);

	# get batch effects or trait vs normalization factor associations
	batch.effects <- NanoStringNorm:::get.batch.effects(x, anno, log, traits, sample.summary.stats);

	# get trait summary stats
	if (check.traits == 1) {
		trait.summary.stats <- NanoStringNorm:::get.trait.summary.stats(x, anno, log, traits);

		# add the trait summary statistics to the gene summary stats
		gene.summary.stats.norm <- cbind(
			gene.summary.stats.norm,
			trait.summary.stats
			);
		}

	# add the normalization details to a list
	x = list(
		normalized.data = data.frame(anno,x),
		raw.data = x.raw,
		normalization.workflow = c(
			CodeCount = CodeCount,
			Background = Background,
			SampleContent = SampleContent,
			otherNorm = otherNorm,
			round = round.values,
			log = log
			),
		sample.summary.stats.norm = sample.summary.stats,
		gene.summary.stats.norm = as.data.frame(gene.summary.stats.norm),
		gene.summary.stats.raw = as.data.frame(gene.summary.stats.raw),
		predicted.concentration = as.data.frame(predicted.concentration),
		batch.effects = batch.effects
		);

	class(x) <- 'NanoStringNorm';
	return(x);

	}
