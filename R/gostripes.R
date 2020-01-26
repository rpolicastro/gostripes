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
		settings = "list"
	),
	prototype(
		sample_sheet = data.frame(),
		settings = list()
	)
)

#' gostripes constructor function
#'
#' Create gostripes object
#'
#' @import methods
#' @import tibble
#'
#' @param sample_sheet Sample sheet data.frame containing 'sample_name', 'replicate_ID", 'R1_read', and optionally 'R2_read'
#' @param cores Number of CPU cores/threads available
#'
#' @rdname gostripes
#'
#' @export

gostripes <- function(sample_sheet, cores = 1) {
	go_obj <- new(
		"gostripes",
		sample_sheet = sample_sheet,
		settings = list("cores" = cores)
	)

	return(go_obj)
}
