
#' Generate Genome Index
#'
#' Generate a STAR genomic index
#'
#' @param go_obj gostripes object
#' @param genome_assembly Genome assembly in fasta file format
#' @param genome_annotation Genome annotation in gtf file format
#' @param outdir Directory to save genome index to
#' @param cores Number of CPU cores/threads to use
#'
#' @rdname genome_index-function
#'
#' @export

genome_index <- function(go_obj, genome_assembly, genome_annotation, outdir, cores = 1) {
	
	## Make sure the output directory exists.
	dir.create(outdir, recursive = TRUE)

	## Store some of the settings.
	go_obj@settings$genome_assembly <- genome_assembly
	go_obj@settings$genome_annotation <- genome_annotation
	go_obj@settings$star_index <- outdir

	## Generate the STAR genome index.
	command <- paste(
		"STAR", "--runMode genomeGenerate",
		"--runThreadN", cores,
		"--genomeDir", outdir,
		"--genomeFastaFiles", genome_assembly,
		"--sjdbGTFfile", genome_annotation
	)
	system(command)

	return(go_obj)
}

#' Align Reads
#'
#' Align reads to genome using STAR
#'
#' @import tibble
#' @importFrom purrr pwalk
#'
#' @param go_obj gostripes object
#' @param outdir Directory to save aligned files to
#' @param cores Number of CPU cores/threads to use
#'
#' @rdname align_reads-function
#'
#' @export

align_reads <- function(go_obj, outdir, cores = 1) {

	## Make sure the output directory exists.
	dir.create(outdir, recursive = TRUE)

	## Store some of the settings.
	go_obj@settings$star_aligned <- outdir

	## Align reads using STAR.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)

		# Get sequencing mode for sample.
		seq_mode <- args$seq_mode

		# Get names of cleaned and dusted fastq files.
		if (seq_mode == "paired") {
			cleaned_R1 <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, "_READ1.fq"))
			cleaned_R2 <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, "_READ2.fq"))
			input_reads <- paste(cleaned_R1, cleaned_R2)
		} else {
			input_reads <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, ".fq"))
		}

		# Align reads using STAR.
		command <- paste(
			"STAR",
			"--runThreadN", cores,
			"--genomeDir", go_obj@settings$star_index,
			"--readFilesIn", input_reads,
			"--outSAMtype BAM SortedByCoordinate",
			"--outFileNamePrefix", file.path(outdir, paste0(args$sample_name, "_"))
		)
		system(command)

		# Index the aligned bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0(args$sample_name, "_Aligned.sortedByCoord.out.bam"))
		)
		system(command)
	})

	return(go_obj)
}
