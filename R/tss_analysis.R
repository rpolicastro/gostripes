
#' Call TSSs
#'
#' @description
#' Call TSSs from BAMs
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
#' @details
#' This function will call TSSs from BAM files.
#' Although during BAM processing there is a tolerance for 3 or less
#' soft-clipped bases, these soft-clipped bases are ignored when
#' calling TSSs.
#'
#' @return gostripes object containing a list of TSSs as GRanges
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   call_TSSs
#'
#' @rdname call_TSSs-function
#'
#' @export

call_TSSs <- function(go_obj) {

	## Check validity of inputs.
	if (!is(go_obj, "gostripes")) stop("go_obj must be gostripes object")

	## Print out some information on TSS calling.
	message("\n## Call TSSs\n")
	
	## Call TSSs for each sample in the sample sheet.
	called_TSSs <- pmap(go_obj@sample_sheet, function(...) {
		args <- list(...)
		message(
			"...Processing ", args$sample_name, "\n",
			"......Calling TSSs"
		)

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
		
		message("......Finished calling TSSs for ", args$sample_name)
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
#' Export TSSs as BEDGRAPHs.
#' The resulting BEDGRAPHs are split into positive and negative stranded TSSs.
#'
#' @importFrom GenomicRanges GRanges strand
#' @importFrom rtracklayer export
#' @importFrom purrr iwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory
#'
#' @return gostripes object and exported TSS BEDGRAPHs
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   call_TSSs %>%
#'   export_TSSs("./scratch/TSSs")
#'
#' @rdname export_TSSs-function
#'
#' @export

export_TSSs <- function(go_obj, outdir) {

	## Check validity of inputs.
	if(!is(go_obj, "gostripes")) stop("go_obj should be a gostripes object")
	if(!is(outdir, "character")) stop("outdir should be a character string")

	## Make sure output directory exists.
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Print out some information on TSS BEDGRAPH export.
	message(
		"\n## TSS BEDGRAPH Export\n",
		"##\n",
		"## Output Directory: ", outdir, "\n"
	)

	## Export TSSs as bedgraphs split by positive and negative strands.
	iwalk(go_obj@TSSs, function(TSSs, sample_name) {
		message(
			"...Processing ", sample_name, "\n",
			"......Calling TSSs"
		)

		# Split TSSs into positive and negative strands.
		pos_TSSs <- TSSs[strand(TSSs) == "+"]
		neg_TSSs <- TSSs[strand(TSSs) == "-"]

		# Create names of bedgraph files.
		pos_bedgraph <- file.path(outdir, paste0("pos_", sample_name, ".bedgraph"))
		neg_bedgraph <- file.path(outdir, paste0("neg_", sample_name, ".bedgraph"))

		# Export bedgraphs.
		export(pos_TSSs, pos_bedgraph, "bedgraph")
		export(neg_TSSs, neg_bedgraph, "bedgraph")

		message("......Finished calling TSSs for ", sample_name)
	})

	return(go_obj)
}

#' Call TSRs
#'
#' @description
#' Basic TSS clustering into TSRs based on naive global read threshold.
#'
#' @import S4Vectors
#' @importFrom GenomicRanges GRanges score reduce
#' @importFrom purrr imap
#'
#' @param go_obj gostripes object
#' @param threshold TSSs with read count below threshold will be discarded
#' @param clust_dist TSSs within this number of base pairs will be clustered
#'
#' @details
#' This function will first remove TSSs with an aggregate read number below 'threshold'.
#' Surviving TSSs will then be clustered if they are within 'clust_dist' bases of eachother.
#'
#' @return gostripes object with list of TSR GRanges
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   call_TSSs %>%
#'   call_TSRs(threshold = 3, clust_dist = 25)
#'
#' @rdname call_TSRs-function
#'
#' @export

call_TSRs <- function(go_obj, threshold, clust_dist) {

	## Check validity of inputs.
	if (!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if (!is(threshold, "numeric")) stop("threshold must be a positive integer")
	if (threshold < 0 | !threshold %% 1 == 0) stop("threshold must be a positive integer")

	## Print out some information on thresholding.
	message(
		"\n## Calling TSRs/cTSSs\n##\n",
		"## Threshold: >= ", threshold, "\n",
		"## Cluster Distance: ", clust_dist, " bases\n"
	)
	
	## Naive thresholding to discover TSS clusters (TSRs or cTSSs).
	TSRs <- imap(go_obj@TSSs, function(TSSs, sample_name) {
		message(
			"...Processing ", sample_name, "\n",
			"......Clustering TSSs into TSRs/cTSSs"
		)

		# Filter TSSs below threshold and merge TSSs within clust_dist.
		filtered_TSSs <- TSSs[score(TSSs) >= threshold]
		clustered <- GenomicRanges::reduce(filtered_TSSs, with.revmap = TRUE, min.gapwidth = clust_dist + 1)

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

		message("......Finished TSR/cTSS calling for ", sample_name)
		return(clustered)
	})

	go_obj@TSRs <- TSRs
	return(go_obj)
}

#' Export TSRs
#'
#' @description
#' Export TSRs as BED files.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer export
#' @importFrom purrr iwalk
#'
#' @param go_obj gostripes object
#' @param outdir Output directory
#'
#' @return gostripes object and TSR BEDs
#'
#' @examples
#' R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
#' R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
#' rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
#' assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
#' annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
#'
#' sample_sheet <- tibble::tibble(
#'   "sample_name" = "stripeseq", "replicate_ID" = 1,
#'   "R1_read" = R1_fastq, "R2_read" = R2_fastq
#' )
#'
#' go_object <- gostripes(sample_sheet) %>%
#'   process_reads("./scratch/cleaned_fastq", rRNA) %>%
#'   fastq_quality("./scratch/fastqc_reports") %>%
#'   genome_index(assembly, annotation, "./scratch/genome_index") %>%
#'   align_reads("./scratch/aligned") %>%
#'   process_bams("./scratch/cleaned_bams") %>%
#'   call_TSSs %>%
#'   call_TSRs(threshold = 3, clust_dist = 25) %>%
#'   export_TSRs("./scratch/TSRs")
#'
#' @rdname export_TSRs-function
#'
#' @export

export_TSRs <- function(go_obj, outdir) {

	## Check validity of inputs.
	if (!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if (!is(outdir, "character")) stop("outdir must be a character string")
	
	## Make sure output directory exists.
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Print out some information on TSR/cTSS export.
	message(
		"\n## TSR/cTSS Export\n##\n",
		"## Output Directory: ", outdir, "\n"
	)

	## Export each TSR as a bed file.
	iwalk(go_obj@TSRs, function(TSRs, sample_name) {		
		message(
			"...Processing ", sample_name, "\n",
			"......Saving TSRs/TSSs to BED"
		)

		# Create new bed file name.
		TSR_bed <- file.path(outdir, paste0(sample_name, ".bed"))

		# Export TSR bed file.
		export(TSRs, TSR_bed, "bed")
		message("......Finished saving TSRs/cTSSs to BED for ", sample_name)
	})

	return(go_obj)
}
