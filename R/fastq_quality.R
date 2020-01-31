
#' Fastq Quality
#'
#' Generate quality reports for FASTQ files
#'
#' @import tibble
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for FastQC reports
#' @param fastq_type Either 'raw', 'processed', or 'both'
#' @param cores Number of CPU cores
#'
#' @rdname fastq_quality-function
#'
#' @export

fastq_quality <- function(go_obj, outdir, fastq_type = "both", cores = 1) {
	
	## Get samples to analyze.
	if (fastq_type %in% c("both", "raw")) {
		raw_R1s <- pull(go_obj@sample_sheet, "R1_read")
		raw_R2s <- go_obj@sample_sheet %>%
			pull("R2_read")
	}
}
