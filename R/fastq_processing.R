
#' Trim R1 Read
#'
#' Check for proper structure of R1 read and then trim it
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet subseq writeXStringSet
#' @importFrom stringr str_which
#' @importFrom purrr pwalk
#'
#' @param go_obj gostripes object
#' @param outdir output directory for trimmed reads
#'
#' @rdname trim_reads-function
#'
#' @export

trim_reads <- function(go_obj, outdir){
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)

		R1_data <- readDNAStringSet(args$R1_read, format = "fastq")
		R1_keep_index <- str_which(R1_data, "^[ATGCN]{8}TATAGGG")
		R1_keep <- R1_data[R1_keep_index] %>%
			subseq(start = 16)

		R2_data <- readDNAStringSet(args$R2_read, format = "fastq")
		R2_keep <- R2_data[R1_keep_index]

		writeXStringSet(R1_keep, file.path(outdir, basename(args$R1_read)), format = "fastq")
		writeXStringSet(R2_keep, file.path(outdir, basename(args$R2_read)), format = "fastq")
	})
}
