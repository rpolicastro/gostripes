#' gostripes class
#'
#' Class container appropriate slots for gostripes
#'
#' @slot sample_sheet Sample sheet
#' @slot settings Settings used throughout gostripes
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
#' @importFrom dplyr mutate
#'
#' @param sample_sheet Sample sheet data.frame containing 'sample_name', 'replicate_ID", 'R1_read', and 'R2_read'
#' @param cores Number of CPU cores/threads available
#'
#' @rdname gostripes
#'
#' @export

gostripes <- function(sample_sheet, cores = 1) {

	## Check whether each sample is paired on unpaired.
	sample_sheet <- sample_sheet %>%
		mutate(seq_mode = ifelse(
			is.na(R2_read) | R2_read %in% c("", " "),
			"unpaired", "paired"
		))

	## Create gostripes object.
	go_obj <- new(
		"gostripes",
		sample_sheet = sample_sheet,
		settings = list("cores" = cores)
	)

	return(go_obj)
}
