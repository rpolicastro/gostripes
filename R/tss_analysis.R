
#' Call TSSs
#'
#' Get TSSs from bams
#'
#' @import tibble
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments
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

		# Check whether paired or sing-end.
		seq_mode <- args$seq_mode

		# Bam file to analyze.
		final_bam <- file.path(go_obj@settings$bam_dir, paste0("soft_", args$sample_name, ".bam"))

		# Read paired-end bam into memory and grab TSSs.
		if (seq_mode == "paired") {
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
				select(-start.first, -end.first)
		# Read single-end bam into memory and grab TSSs.
		} else {
			TSSs <- final_bam %>%
				readGAlignments(use.name = TRUE) %>%
				as.data.frame %>%
				as_tibble(.name_repair = "unique", rownames = "qname") %>%
				select(seqnames, strand, start, end) %>%
				mutate(
					"start" = ifelse(strand == "+", start, end),
					"end" = start
				)
		}

		# Aggregate and sum overlapping TSSs.
		TSSs <- TSSs %>%
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

#' Export TSSs
#'
#' Export TSSs as bedgraphs
#'
#' @importFrom GenomicRanges GRanges strand
#' @importFrom rtracklayer export
#' @importFrom purrr iwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory
#'
#' @rdname export_TSSs-function
#'
#' @export

export_TSSs <- function(go_obj, outdir) {

	## Make sure output directory exists.
	if (!dir.exists(outdir)) {
		dir.create(outdir, recursive = TRUE)
	}

	## Export TSSs as bedgraphs split by positive and negative strands.
	iwalk(go_obj@TSSs, function(TSSs, sample_name) {

		# Split TSSs into positive and negative strands.
		pos_TSSs <- TSSs[strand(TSSs) == "+"]
		neg_TSSs <- TSSs[strand(TSSs) == "-"]

		# Create names of bedgraph files.
		pos_bedgraph <- file.path(outdir, paste0("pos_", sample_name, ".bedgraph"))
		neg_bedgraph <- file.path(outdir, paste0("neg_", sample_name, ".bedgraph"))

		# Export bedgraphs.
		export(pos_TSSs, pos_bedgraph, "bedgraph")
		export(neg_TSSs, neg_bedgraph, "bedgraph")
	})
}

#' Call TSRs
#'
#' Basic TSS clustering into TSRs based on naive global read threshold
#'
#' @import S4Vectors
#' @importFrom GenomicRanges GRanges score reduce
#' @importFrom purrr map
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

		# Get the aggregated score and unique TSSs of clustered TSSs.
		cluster_info <- aggregate(
			filtered_TSSs,
			mcols(clustered)$revmap,
			score = sum(score),
			n_unique = lengths(score)
		)

		# Add aggregated scores and unique TSS numbers back to clustered TSSs.
		clustered$score <- cluster_info$score
		clustered$n_unique <- cluster_info$n_unique
		clustered$revmap <- NULL

		return(clustered)
	})

	go_obj@TSRs <- TSRs
	return(go_obj)
}

#' Export TSRs
#'
#' Export TSRs as bed files
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer export
#' @importFrom purrr iwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory
#'
#' @rdname export_TSRs-function
#'
#' @export

export_TSRs <- function(go_obj, outdir) {
	
	## Make sure output directory exists.
	if (!dir.exists(outdir)) {
		dir.create(outdir, recursive = TRUE)
	}

	## Export each TSR as a bed file.
	iwalk(go_obj@TSRs, function(TSRs, sample_name) {
		
		# Create new bed file name.
		TSR_bed <- file.path(outdir, paste0(sample_name, ".bed"))

		# Export TSR bed file.
		export(TSRs, TSR_bed, "bed")
	})
}
