
#' Process BAM Files
#'
#' @description
#' Quality control steps for BAM files.
#' PCR duplicates are first removed along with non-primary alignments,
#' and then reads missing mate pairs (for paired-end data).
#' R1 reads with more than 3 soft-clipped bases on the 5' end of the read
#' are then discarded.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom purrr pwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for cleaned reads
#' @param cores Number of CPU cores/threads to use
#'
#' @details
#' For paired-end data, Samtools is used both to mark and remove duplicates,
#' as well as remove non-primary alignments and reads missing mate pairs.
#' This produces BAMs in the output directory with the names 'deduped_*'.
#'
#' For single-end data, samtools is used to remove non-primary reads,
#' producing BAMs with the name 'filtered_*'.
#' Then, UMI-tools is used to remove PCR duplicates via the previously stashed UMIs.
#' The file name for these BAMs is 'deduped_*'.
#'
#' After the above processing steps, R1 reads with more than 3 soft-clipped 5' bases
#' are removed. Some tolerance is allowed for soft-clipped bases because it is common for
#' there to be at least one non-templated C adjacent to the true TSS perportedly due to
#' the cap acting as a template itself.
#'
#' @return gostripes object and processed BAM files
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
#'   process_bams("./scratch/cleaned_bams")
#'
#' @rdname process_bams-function
#'
#' @export

process_bams <- function(go_obj, outdir, cores = 1) {

	## Check validity of inputs.
	if (!is(go_obj, "gostripes")) stop("go_obj should be a gostripes object")
	if (!is(outdir, "character")) stop("outdir should be a character string")
	if (!is(cores, "numeric")) stop("cores should be a positive integer")
	if (cores < 1 | !cores %% 1 == 0) stop("cores should be a positive integer")

	## Ensure output directory exists.
	dir.create(outdir, recursive = TRUE)

	## Saving some settings.
	go_obj@settings$bam_dir <- outdir

	# Print out some information about the bam processing.
	message(
		"\n## Bam Processing\n",
		"##\n",
		"## Output Directory: ", outdir, "\n",
		"## Cores: ", cores, "\n"
	)

	## Processing the BAM files.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)
		message("...Processing ", args$sample_name)

		# Remove PCR duplicates, reads missing mate, and non-primary alignments.
		bam_file <- file.path(go_obj@settings$star_aligned, paste0(args$sample_name, "_Aligned.sortedByCoord.out.bam"))

		message("......Removing PCR duplicates, non-primary alignments, and for paired-end data reads missing mates")
		filter_flags(bam_file, args$sample_name, outdir, args$seq_mode, cores)

		# Remove reads where the R1 has more than 3 soft-clipped bases on the 5' end.
		deduped_bam <- file.path(outdir, paste0("deduped_", args$sample_name, ".bam"))

		message("......Removing fragments with too many soft-clipped bases near TSS")
		filter_soft(deduped_bam, args$sample_name, outdir, args$seq_mode)

		message("......Done processing ", args$sample_name, "!")
	})

	return(go_obj)
}

#' Filter Flags
#'
#' Remove PCR duplicates and other undesirable reads
#'
#' @param bam_file Aligned bam file to process
#' @param sample_name Name of sample
#' @param outdir Output directory
#' @param seq_mode Whether the run was paired-end or single-end
#' @param cores Number of CPU cores available
#'
#' @rdname filter_flag-function
#'
#' @export

filter_flags <- function(bam_file, sample_name, outdir, seq_mode, cores) {

	## Remove PCR duplicates and undesirable flags.
	if (seq_mode == "paired") {

		# For paired-end data use samtools to remove both duplicates and undesirable flags.
		command <- paste(
			"samtools sort -n -@", cores, bam_file, "|",
			"samtools fixmate -m - - |",
			"samtools sort -@", cores, "- |",
			"samtools markdup - - |",
			"samtools view -F 3852 -f 3 -O BAM -@", cores,
			"-o", file.path(outdir, paste0("deduped_", sample_name, ".bam"))
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("deduped_", sample_name, ".bam"))
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

	} else {

		# For single-end data use samtools to just remove undesirable flags.
		command <- paste(
			"samtools view -F 3844 -O BAM -@", cores,
			"-o", file.path(outdir, paste0("filtered_", sample_name, ".bam")),
			bam_file
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("filtered_", sample_name, ".bam"))
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

		# Use UMI-tools to remove PCR duplicates.
		command <- paste(
			"umi_tools dedup",
			"-I", file.path(outdir, paste0("filtered_", sample_name, ".bam")),
			"-S", file.path(outdir, paste0("deduped_", sample_name, ".bam")),
			paste0("--output-stats=", file.path(outdir, paste0("dedup_", sample_name)))
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("deduped_", sample_name, ".bam"))
		)
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
	}

	
}

#' Filter Soft-clipped
#'
#' Filter reads with too many soft-clipped bases
#'
#' @import tibble
#' @importFrom dplyr select mutate filter rename
#' @importFrom stringr str_extract
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
#'
#' @param flagfiltered_bam Bam file with PCR duplicates and undesirable flags filtered out
#' @param sample_name Name of the sample
#' @param outdir Output directory
#' @param seq_mode Whether the bam is paired-end or single-end
#'
#' @rdname filter_soft-function
#'
#' @export

filter_soft <- function(deduped_bam, sample_name, outdir, seq_mode) {

	## Load R1 reads from bam.
	if (seq_mode == "paired") {
		# Read in pairs if paired-end.
		bam_data <- deduped_bam %>%
			readGAlignmentPairs(use.names = TRUE) %>%
			as.data.frame %>%
			as_tibble(.name_repair = "unique", rownames = "qname") %>%
			select(qname, strand.first, cigar.first) %>%
			rename(strand = strand.first, cigar = cigar.first)
	} else {
		# Read in single end if not paired.
		bam_data <- deduped_bam %>%
			readGAlignments(use.names = TRUE) %>%
			as.data.frame %>%
			as_tibble(.name_repair = "unique", rownames = "qname") %>%
			select(qname, strand, cigar)
	}

	## Retrieve R1 reads with more than 3 soft-clipped 5' bases.
	R1_keep <- bam_data %>%
		mutate(
			cigar_soft = ifelse(
				strand == "-",
				str_extract(cigar, "\\d+S$"),
				str_extract(cigar, "^\\d+S")
			),
			n_soft = str_extract(cigar_soft, "^\\d+") %>% as.numeric
		) %>%
		filter(!is.na(n_soft) & n_soft > 3) %>%
		select(qname)

	## Write read names to exclude to text file.
	write.table(
		R1_keep, file.path(outdir, paste0("soft_", sample_name, ".txt")),
		col.names = FALSE, row.names = FALSE, quote = FALSE
	)

	## Construct picard tools command to remove the reads marked with too many soft clipped bases.
	command = paste(
		"picard FilterSamReads",
		paste0("I=", deduped_bam),
		paste0("O=", file.path(outdir, paste0("soft_", sample_name, ".bam"))),
		paste0("READ_LIST_FILE=", file.path(outdir, paste0("soft_", sample_name, ".txt"))),
		"FILTER=excludeReadList"
	)
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

	## Index the bams.
	command <- paste(
		"samtools index",
		file.path(outdir, paste0("soft_", sample_name, ".bam"))
	)
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
}
