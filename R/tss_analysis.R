
#' Call TSSs
#'
#' Get TSSs from bams
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom GenomicAlignments readGAlignmentPairs GAlignmentPairs
#' @importFrom dplyr select rename mutate count pull
#' @importFrom purrr pmap
#' @importFrom magrittr %>% set_names
#'
#' @param go_obj gostripes object
#'
#' @rdname call_TSSs-function
#'
#' @export

call_TSSs <- function(go_obj) {
	
	## Call TSSs for each sample in the sample sheet.
	called_TSSs <- pmap(go_obj@sample_sheet, function(...) {
		args <- list(...)

		# Bam file to analyze.
		final_bam <- file.path(go_obj@settings$bam_dir, paste0("final_", args$sample_name, ".bam"))

		# Read bam into memory and call TSSs.
		TSSs <- final_bam %>%
			readGAlignmentPairs(use.names = TRUE) %>%
			as.data.frame %>%
			as_tibble(.name_repair = "unique", rownames = "qname") %>%
			select(seqnames.first, strand.first, start.first, end.first) %>%
			mutate(
				"start" = ifelse(strand == "+", start.first, end.first),
				"end" = start
			) %>%
			rename("seqnames" = seqnames.first, "strand" = strand.first) %>%
			select(-start.first, -end.first) %>%
			count(seqnames, start, end, strand, name = "score") %>%
			makeGRangesFromDataFrame(keep.extra.columns = TRUE)

		return(TSSs)
	})

	## Add names to the list of TSS Granges.
	called_TSSs <- set_names(called_TSSs, pull(go_obj@sample_sheet, "sample_name"))

	## Add TSSs back to gostripes object.
	go_obj@TSSs <- called_TSSs

	return(go_obj)
}
