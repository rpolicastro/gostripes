
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

		filter_flags(bam_file, args$sample_name, outdir, cores)

		# Remove read pairs where the R1 has more than 3 soft-clipped bases on the 5' end.
		flag_filtered_bam <- file.path(outdir, paste0("flagfiltered_", args$sample_name, ".bam"))

		filter_soft(flag_filtered_bam, args$sample_name, outdir)
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
#' @param cores Number of CPU cores available
#'
#' @rdname filter_flag-function
#'
#' @export

filter_flags <- function(bam_file, sample_name, outdir, cores) {

	## Make samtools command to mark and remove duplicates, non-primary reads, and reads missing mate pairs.
	command <- paste(
		"samtools sort -n -@", cores, bam_file, "|",
		"samtools fixmate -m - - |",
		"samtools sort -@", cores, "- |",
		"samtools markdup - - |",
		"samtools view -F 3852 -f 3 -O BAM -@", cores,
		"-o", file.path(outdir, paste0("flagfiltered_", sample_name, ".bam"))
	)

	## Run samtools command.
	system(command)
	
}

#' Filter Soft-clipped
#'
#' Filter reads with too many soft-clipped bases
#'
#' @import tibble
#' @importFrom dplyr select mutate filter
#' @importFrom stringr str_extract
#' @importFrom GenomicAlignments readGAlignmentPairs 
#'
#' @param flagfiltered_bam Bam file with PCR duplicates and undesirable flags filtered out
#' @param sample_name Name of the sample
#' @param outdir Output directory
#'
#'
#' @rdname filter_soft-function
#'
#' @export

filter_soft <- function(flagfiltered_bam, sample_name, outdir) {

	## Load R1 reads from bam.
	bam_pairs <- flagfiltered_bam %>%
		readGAlignmentPairs(use.names = TRUE) %>%
		as.data.frame %>%
		as_tibble(.name_repair = "unique", rownames = "qname") %>%
		select(qname, strand.first, cigar.first)

	## Retrieve R1 reads with more than 3 soft-clipped 5' bases.
	R1_keep <- bam_pairs %>%
		mutate(
			cigar_soft = ifelse(
				strand.first == "-",
				str_extract(cigar.first, "\\d+S$"),
				str_extract(cigar.first, "^\\d+S")
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
		paste0("I=", flagfiltered_bam),
		paste0("O=", file.path(outdir, paste0("final_", sample_name, ".bam"))),
		paste0("READ_LIST_FILE=", file.path(outdir, paste0("softfiltered_", sample_name, ".txt"))),
		"FILTER=excludeReadList"
	)

	## Run the picard tools command.
	system(command)
	
}
