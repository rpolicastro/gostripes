
#' Count Features
#'
#' Count reads associated with annotated features
#'
#' @import tibble
#' @importFrom Rsubread featureCounts
#' @importFrom purrr map
#' @importFrom dplyr pull
#'
#' @param go_obj gostripes object
#' @param genome_annotation Genome annotation in GTF file format
#' @param cores Number of CPU cores available
#'
#' @rdname count_features-function
#'
#' @export

count_features <- function(go_obj, genome_annotation, cores = 1) {

	## Separate paired-end and single-end reads.
	seq_status <- go_obj@sample_sheet %>%
		split(.$seq_mode) %>%
		map(
			~ pull(., "sample_name") %>%
			file.path(go_obj@settings$bam_dir, paste0("soft_", ., ".bam"))
		)

	## Build featureCounts command.
	counts <- imap(seq_status, function(bams, seq_mode) {

		# Count paired-end features.
		if (seq_mode == "paired") {
			feature_counts <- featureCounts(
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
			)

		# Count single-end features.
		} else {
			feature_counts <- featureCounts(
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
			)
		}

		return(feature_counts)
	})

	## Merge counts.
}
