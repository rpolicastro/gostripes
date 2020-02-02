
#' Count Features
#'
#' @description
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
#' @details
#' Genome annotations can be found in repositories such as NCBI, UCSC, and ensembl.
#' The 'genome_annotation' file should be the same GTF used in read alignment for consistency.
#'
#' This function uses the featureCounts function from Rsubread to summarize counts to annotated features.
#' First, sequenced fragments are assigned to the nearest exon if there is at least 10 overlapping bases.
#' If the fragment overlaps more than one feature, it is asigned to the feature with the largest overlap.
#' Finally, counts for all exons from the same gene are aggregated into a sum score for that gene.
#'
#' @return gostripes object with feature counts matrix
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   count_features(annotation)
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
#' @description
#' Export feature counts as a table
#'
#' @import tibble
#' @importFrom dplyr pull
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for table
#'
#' @return gostripes object and tab separated table of feature counts.
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   count_features(annotation) %>%
#'   export_counts("./scratch/counts")
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
		"## Output Directory: ", outdir, "\n"
	)

	## Export the counts to a table.
	message("...Exporting feature counts table")
	write.table(
		go_obj@feature_counts, file.path(outdir, "feature_counts.tsv"),
		col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
	)
	message("...Finished exporting feature counts table")

	return(go_obj)
}
