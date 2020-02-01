
#' Count Features
#'
#' Count reads associated with annotated features
#'
#' @import tibble
#' @importFrom Rsubread featureCounts
#' @importFrom purrr map imap reduce
#' @importFrom dplyr pull left_join
#' @importFrom stringr str_replace
#'
#' @param go_obj gostripes object
#' @param genome_annotation Genome annotation in GTF file format
#' @param cores Number of CPU cores available
#'
#' @rdname count_features-function
#'
#' @export

count_features <- function(go_obj, genome_annotation, cores = 1) {

	## Print out some information on feature counting.
	message(
		"\n## Feature Counting\n",
		"##\n",
		"## Annotation: ", genome_annotation, "\n",
		"## Cores: ", cores, "\n",
		"##\n",
		"## Feature Count Settings:\n",
		"## - Assign fragments to exon \n",
		"## - Summarize fragment counts by gene \n",
		"## - Fragments must overlap feature at least 10 bases \n",
		"## - Fragments are assigned to feature with largest overlap \n",
		"##\n",
		"## Samples to Analyze:\n",
		sprintf("## - %s\n", pull(go_obj@sample_sheet, "sample_name")), "\n\n",
		"...Started counting features"
	)
		

	## Separate paired-end and single-end reads.
	seq_status <- go_obj@sample_sheet %>%
		split(.$seq_mode) %>%
		map(function(x) {
			samp_names <- pull(x, "sample_name")
			samp_names <- file.path(go_obj@settings$bam_dir, paste0("soft_", samp_names, ".bam"))
			return(samp_names)
		})

	## Build featureCounts command.
	counts <- imap(seq_status, function(bams, seq_mode) {

		# Count paired-end features.
		if (seq_mode == "paired") {
			capture.output(feature_counts <- featureCounts(
				files = bams,
				annot.ext = genome_annotation,
				isGTFAnnotationFile = TRUE,
				GTF.featureType = "exon",
				GTF.attrType = "gene_id",
				useMetaFeatures = TRUE,
				allowMultiOverlap = FALSE,
				minOverlap = 10,
				largestOverlap = TRUE,
				strandSpecific = 1,
				isPairedEnd = TRUE,
				nthreads = cores
			))

		# Count single-end features.
		} else {
			capture.output(feature_counts <- featureCounts(
				files = bams,
				annot.ext = genome_annotation,
				isGTFAnnotationFile = TRUE,
				GTF.featureType = "exon",
				GTF.attrType = "gene_id",
				useMetaFeatures = TRUE,
				allowMultiOverlap = FALSE,
				minOverlap = 10,
				largestOverlap = TRUE,
				strandSpecific = 1,
				isPairedEnd = FALSE,
				nthreads = cores,
				readExtension3 = 200
			))
		}

		# Extract feature counts and remove .bam from sample names.
		feature_counts <- feature_counts$counts %>%
			as_tibble(.name_repair = "unique", rownames = "gene_id")
		colnames(feature_counts) <- str_replace(colnames(feature_counts), "\\.bam$", "")

		return(feature_counts)
	})

	## Merge counts.
	counts <- reduce(counts, left_join, by = "gene_id")
	message("...Finished counting features!")

	## Add counts back to gostripes object.
	go_obj@feature_counts <- counts
	return(go_obj)
}

#' Export Feature Counts
#'
#' Export feature counts as a table
#'
#' @import tibble
#' @importFrom dplyr pull
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for table
#'
#' @rdname export_counts-function
#'
#' @export

export_counts <- function(go_obj, outdir) {

	## Check validity of inputs.
	if(!is(go_object, "gostripes")) stop("go_obj should be a gostripes object")
	if(!is(outdir, "character")) stop("outdir should be a character string")

	## Ensure output directory exists.
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Print out some information.
	message(
		"\n## Exporting Feature Counts\n",
		"##\n",
		"## Output Directory: ", outdir, "\n",
		"##\n",
		"## Samples with Feature Counts:\n",
		sprintf("## - %s\n", pull(go_obj@sample_sheet, "sample_name"))
	)

	## Export the counts to a table.
	message("...Exporting feature counts table")
	write.table(
		go_obj@feature_counts, file.path(outdir, "feature_counts.tsv"),
		col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
	)
	message("...Finished exporting feature counts table")
}
