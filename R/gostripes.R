#' gostripes class
#'
#' Class container appropriate slots for gostripes
#'
#' @slot sample_sheet Sample sheet
#' @slot settings Settings used throughout gostripes
#' @slot TSSs List of GRanges containing transcription start sites
#' @slot TSRs List of GRagnes containing transcription start regions (TSRs) or clustered transcription start sites (cTSSs)
#' @slot feature_counts data.frame containing the RNA-seq like feature counting of reads
#'
#' @rdname gostripes

setClass(
	"gostripes",
	representation(
		sample_sheet = "data.frame",
		settings = "list",
		TSSs = "list",
		TSRs = "list",
		feature_counts = "data.frame"
	),
	prototype(
		sample_sheet = data.frame(),
		settings = list(),
		TSSs = list(),
		TSRs = list(),
		feature_counts = data.frame()
	)
)

#' gostripes constructor function
#'
#' Create gostripes object
#'
#' @import methods
#' @import tibble
#' @importFrom dplyr mutate pull
#' @importFrom magrittr %>%
#'
#' @param sample_sheet Sample sheet data.frame containing 'sample_name', 'replicate_ID", 'R1_read', and 'R2_read'
#' @param cores Number of CPU cores/threads available
#'
#' @rdname gostripes
#'
#' @export

gostripes <- function(sample_sheet, cores = 1) {

	## Check for proper inputs to gostripes function.

	# Check whether sample_sheet is in proper format.
	if (!is(sample_sheet, "data.frame")) stop("The sample sheet must be a data frame or tibble.")
	if (nrow(sample_sheet) == 0) stop("The sample sheet contains no entries.")
	if (!all(c("sample_name", "replicate_ID", "R1_read", "R2_read") %in% colnames(sample_sheet))) {
		stop("The sample sheet must have columns 'sample_name', 'replicate_ID', 'R1_read', and 'R2_read'.")
	}

	# Check other arguments.
	if (!is(cores, "numeric")) stop("Cores must be a positive integer")
	if (!cores %% 1 == 0 | cores < 1) stop("Cores must be a positive integer")

	## Check whether each sample is paired on unpaired.
	sample_sheet <- sample_sheet %>%
		mutate(seq_mode = ifelse(
			is.na(R2_read) | R2_read %in% c("", " "),
			"unpaired", "paired"
		))

	## Print out some information on the sample sheet.
	message(
		"## gostripesR v0.2.0\n##\n",
		"## Sample sheet contains ", nrow(sample_sheet), " sample(s)\n",
		sprintf("## - %s\n", pull(sample_sheet, "sample_name")),
		"##\n",
		"## Available cores set to ", cores
	)

	## Create gostripes object.
	go_obj <- new(
		"gostripes",
		sample_sheet = sample_sheet,
		settings = list("cores" = cores)
	)

	return(go_obj)
}
