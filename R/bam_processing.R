
#' Process BAM Files
#'
#' Quality control steps for BAM files.
#' Includes PCR duplicate removal, ensuring both pairs mapped, and keeping only primary alignments.
#'
#' @import tibble
#' @importFrom magrittr %>%
#' @importFrom purrr pwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for cleaned reads
#' @param cores Number of CPU cores/threads to use
#'
#' @rdname process_bams-function
#'
#' @export

process_bams <- function(go_obj, outdir, cores = 1) {

	## Ensure output directory exists.
	dir.create(outdir, recursive = TRUE)

	## Saving some settings.
	go_obj@settings$bam_dir <- outdir

	## Processing the BAM files.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)

		# Remove PCR duplicates, reads missing mate, and non-primary alignments.
		bam_file <- file.path(go_obj@settings$star_aligned, paste0(args$sample_name, "_Aligned.sortedByCoord.out.bam"))

		filter_flags(bam_file, args$sample_name, outdir, args$seq_mode, cores)

		# Remove reads where the R1 has more than 3 soft-clipped bases on the 5' end.
		deduped_bam <- file.path(outdir, paste0("deduped_", args$sample_name, ".bam"))

		filter_soft(deduped_bam, args$sample_name, outdir, args$seq_mode)
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
		system(command)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("deduped_", sample_name, ".bam"))
		)
		system(command)

	} else {

		# For single-end data use samtools to just remove undesirable flags.
		command <- paste(
			"samtools view -F 3844 -O BAM -@", cores,
			"-o", file.path(outdir, paste0("filtered_", sample_name, ".bam")),
			bam_file
		)
		system(command)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("filtered_", sample_name, ".bam"))
		)
		system(command)

		# Use UMI-tools to remove PCR duplicates.
		command <- paste(
			"umi_tools dedup",
			"-I", file.path(outdir, paste0("filtered_", sample_name, ".bam")),
			"-S", file.path(outdir, paste0("deduped_", sample_name, ".bam")),
			paste0("--output-stats=", file.path(outdir, paste0("dedup_", sample_name)))
		)
		system(command)

		# Index the bams.
		command <- paste(
			"samtools index",
			file.path(outdir, paste0("deduped_", sample_name, ".bam"))
		)
		system(command)
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
		R1_keep, file.path(outdir, paste0("softfiltered_", sample_name, ".txt")),
		col.names = FALSE, row.names = FALSE, quote = FALSE
	)

	## Construct picard tools command to remove the reads marked with too many soft clipped bases.
	command = paste(
		"picard FilterSamReads",
		paste0("I=", deduped_bam),
		paste0("O=", file.path(outdir, paste0("softfiltered_", sample_name, ".bam"))),
		paste0("READ_LIST_FILE=", file.path(outdir, paste0("softfiltered_", sample_name, ".txt"))),
		"FILTER=excludeReadList"
	)
	system(command)

	## Index the bams.
	command <- paste(
		"samtools index",
		file.path(outdir, paste0("softfiltered_", sample_name, ".bam"))
	)
	system(command)
}
