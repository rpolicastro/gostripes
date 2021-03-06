
#' Fastq Quality
#'
#' @description
#' Generate quality control reports for FASTQ files using FastQC
#'
#' @import tibble
#' @importFrom dplyr filter pull mutate
#' @importFrom purrr pmap reduce
#'
#' @param go_obj gostripes object
#' @param outdir Output directory for FastQC reports
#' @param fastq_type Either 'raw', 'processed', or 'both'
#' @param cores Number of CPU cores
#'
#' @details
#' The fastq_type specifies the FASTQ files to quality control.
#' 'raw' will analyze the FASTQ files before the initial processing
#' and quality control steps. 'processed' will look at the FASTQ files
#' after the quality control and processing steps. 'both' will generate
#' reports for the 'raw' and 'processed' FASTQ files.
#'
#' @return FastQC quality reports
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports")
#'
#' @rdname fastq_quality-function
#'
#' @export

fastq_quality <- function(go_obj, outdir, fastq_type = "both", cores = 1) {

	## Check inputs to argument.
	if (!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if (!is(outdir, "character")) stop("outdir must be a character")
	if (!is(fastq_type, "character")) stop("fastq_type must be 'raw', 'processed', or 'both'")
	if (!fastq_type %in% c("both", "raw", "processed")) stop("fastq_type must be 'raw', 'processed', or 'both'")
	if (!is(cores, "numeric")) stop("cores must be a positive integer")
	if (!cores %% 1 == 0 | cores < 1) stop("cores must be a positive integer")

	## Make sure output directory exists.
	if (!dir.exists(outdir)) {
		dir.create(outdir, recursive = TRUE)
	}
	
	## Select raw R1 and also R2 reads if paired-end.
	if (fastq_type %in% c("both", "raw")) {
		raw_R1s <- pull(go_obj@sample_sheet, "R1_read")
		raw_R2s <- go_obj@sample_sheet %>%
			filter(seq_mode == "paired") %>%
			pull("R2_read")

		if (length(raw_R2s) == 0) {
			raw_fastqs <- raw_R1s
		} else {
			raw_fastqs <- c(raw_R1s, raw_R2s)
		}
	}

	## Select processed fastq files if required.
	if (fastq_type %in% c("both", "processed")) {
		processed_fastqs <- pmap(go_obj@sample_sheet, function(...) {
			args <- list(...)			

			# Get sequencing mode.
			seq_mode <- args$seq_mode

			# Prepare R1 read.
			processed_R1 <- ifelse(
				seq_mode == "paired",
				file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, "_READ1.fq")),
				file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, ".fq"))
			)

			# Prepare R2 read if paired end.
			if (seq_mode == "paired") {
				processed_R2 <- file.path(
					go_obj@settings$fastq_outdir,
					paste0("decon_", args$sample_name, "_READ2.fq")
				)
			}

			# Return the list of processed fastq files.
			if (seq_mode == "paired") {
				processed_fastqs <- c(processed_R1, processed_R2)
			} else {
				processed_fastqs <- processed_R1
			}
		
			return(processed_fastqs)
		}) %>%
		reduce(c)
	}

	## Create complete list of FASTQ files to analyze.
	if (fastq_type == "raw") {
		fastqs <- raw_fastqs
	} else if (fastq_type == "processed") {
		fastqs <- processed_fastqs
	} else {
		fastqs <- c(raw_fastqs, processed_fastqs)
	}

	## Print out some information on the FastQC analysis.
	message(
		"\n## FastQC Analysis of Raw and/or Processed Reads\n",
		"##\n",
		"## Output Directory: ", outdir, "\n",
		"## Cores: ", cores, "\n\n",
		"...Started FastQC quality control"
	)

	## Run FastQC quality control on FASTQ files.
	command <- paste("fastqc -t", cores, "-o", outdir, paste0(fastqs, collapse = " "))
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
	message("...Finished FastQC quality control!")

	go_obj@settings$fastqc_outdir <- outdir
	return(go_obj)
}
