
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
				"start" = ifelse(strand.first == "+", start.first, end.first),
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

#' Call TSRs
#'
#' Basic TSS clustering into TSRs based on naive global read threshold
#'
#' @import S4Vectors
#' @importFrom GenomicRanges GRanges score reduce
#'
#' @param go_obj gostripes object
#' @param threshold TSSs with read count below threshold will be discarded
#' @param clust_dist TSSs within this number of base pairs will be clustered
#'
#' @rdname call_TSRs-function
#'
#' @export

call_TSRs <- function(go_obj, threshold, clust_dist) {
	
	## Naive thresholding to discover TSS clusters (TSRs).
	TSRs <- map(go_obj@TSSs, function(TSSs) {

		# Filter TSSs below threshold and merge TSSs within clust_dist.
		filtered_TSSs <- TSSs[score(TSSs) >= threshold]
		clustered <- reduce(filtered_TSSs, with.revmap = TRUE, min.gapwidth = clust_dist + 1)

		# Get the aggregated score of clustered TSSs.
		cluster_scores <- aggregate(filtered_TSSs, mcols(clustered)$revmap, score = sum(score))

		# Add aggregated scores back to clustered TSSs.
		clustered$score <- cluster_scores$score
		clustered$revmap <- NULL

		return(clustered)
	})

	go_obj@TSRs <- TSRs
	return(go_obj)
}
