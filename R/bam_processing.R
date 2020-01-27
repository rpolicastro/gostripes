
#' Process BAM Files
#'
#' Quality control steps for BAM files.
#' Includes PCR duplicate removal, ensuring both pairs mapped, and keeping only primary alignments.
#'
#' @import tibble
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom dplyr pull select rename mutate filter
#' @importFrom stringr str_extract
#' @importFrom magrittr %>%
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for cleaned reads
#' @param cores Number of CPU cores/threads to use
#'
#' @rdname process_bams-function
#'
#' @export

process_bams <- function(go_obj, outdir, cores = 1) {

	## Saving some settings.
	go_obj@settings$bam_dir <- outdir

	## Processing the BAM files.
	pwalk(go_obj@sample_sheet, function(...) {
		# File name of aligned reads.
		bam_file <- file.path(deep_obj@settings$star_aligned, paste0(args$sample_name, "_Aligned.out.bam"))

		# Create samtools command to remove duplicates, non-primary reads, and reads without mate.
		command <- paste(
			"samtools sort -n -@", cores, bam_file, "|",
			"samtools fixmate -m - - |",
			"samtools sort -@", cores, "- |",
			"samtools markdup - - |",
			"samtools view -F 3852 -f 3 -O BAM -@", cores,
			"-o", file_path(outdir, paste0(args$sample_name, ".bam"))
		)
		system(command)

		# Load R1 reads from bam.
		bam_pairs <- file_path(outdir, paste0(args$sample_name, ".bam")) %>%
			readGAlignmentPairs(use.names = TRUE) %>%
			as.data.frame %>%
			as_tibble(.name_repair = "unique", rownames = "qname") %>%
			select(qname, strand.first, cigar.first)

		# Retrieve R1 reads with 3 or less soft-clipped 5' bases.
		R1_keep <- bam_pairs %>%
			mutate(
				cigar_soft = ifelse(
					strand.first == "-",
					str_extract(cigar.first, "\\d+S$"),
					str_extract(cigar.first, "^\\d+S")
				),
				n_soft = str_extract(cigar_soft, "^\\d+") %>% as.numeric		
			) %>%
			filter(is.na(n_soft) | n_soft <= 3) %>%
			pull(qname)
	})

	return(go_obj)
}
