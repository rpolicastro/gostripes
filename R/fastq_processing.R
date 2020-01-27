
#' Process Reads
#'
#' Remove rRNA contamination, low complexity reads, and trim UMI-Spacer-GGG
#'
#' @import tibble
#' @importFrom Biostrings readDNAStringSet DNAStringSet subseq writeXStringSet
#' @importFrom stringr str_which
#' @importFrom purrr pwalk
#' @importFrom magrittr %>%
#'
#' @param go_obj gostripes object
#' @param contamination_fasta fasta file containing contaminants to remove, such as rRNA
#' @param outdir output directory for filtered and trimmed reads
#' @param cores Number of CPU core/threads to use
#'
#' @rdname process_reads-function
#'
#' @export

process_reads <- function(go_obj, outdir, contamination_fasta, cores = 1){
	
	## Add output directory to settings.
	go_obj@settings$fastq_outdir <- outdir

	## Process each pair of reads.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)

		# Load up R1 read and check for reads matching proper R1 structure.
		R1_data <- readDNAStringSet(args$R1_read, format = "fastq")
		R1_keep_index <- str_which(R1_data, "^[ATGCN]{8}TATAGGG")
		R1_keep <- R1_data[R1_keep_index] %>%
			subseq(start = 16)

		# Discard the matching paired end R2 read.
		R2_data <- readDNAStringSet(args$R2_read, format = "fastq")
		R2_keep <- R2_data[R1_keep_index]

		# Write trimmed fastq read to file.
		writeXStringSet(R1_keep, file.path(outdir, basename(args$R1_read)), format = "fastq")
		writeXStringSet(R2_keep, file.path(outdir, basename(args$R2_read)), format = "fastq")

		# TagDust2 to remove contaminants such as rRNA and low complexity reads.
		command <- paste(
			"tagdust",
			"-ref", contamination_fasta,
			"-fe 3", "-t", cores, "-dust 97",
			"-o", file.path(outdir, args$sample_name), "-1 R:N",
			file.path(outdir, basename(args$R1_read)),
			file.path(outdir, basename(args$R2_read))
		)
		system(command)
	})

	return(go_obj)
}
