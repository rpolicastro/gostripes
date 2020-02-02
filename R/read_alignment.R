
#' Generate Genome Index
#'
#' @description
#' Generate a STAR genomic index to be used in read alignment
#'
#' @param go_obj gostripes object
#' @param genome_assembly Genome assembly in FASTA file format
#' @param genome_annotation Genome annotation in GTF file format
#' @param outdir Directory to save genome index to
#' @param cores Number of CPU cores/threads to use
#'
#' @details
#' Genome assemblies and annotations can be found in repositories such as
#' NCBI, UCSC, and ensembl. The 'genome_assembly' should be a properly formatted
#' FASTA file, ideally with repeats soft-masked. The 'genome_annotation' file
#' should be in GTF format, as the STAR aligner does not use the newer GFF3 format.
#'
#' @return A STAR geome index
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
#'   genome_index(assembly, annotation, "./scratch/genome_index")
#'
#' @rdname genome_index-function
#'
#' @export

genome_index <- function(go_obj, genome_assembly, genome_annotation, outdir, cores = 1) {

	## Check validity of inputs.
	if (!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if (!is(genome_assembly, "character")) stop("genome_assembly must be a character string")
	if (!file.exists(genome_assembly)) stop("genome_assembly file does not exist")
	if (!is(genome_annotation, "character")) stop("genome_annotation must be a character string")
	if (!file.exists(genome_annotation)) stop("genome_annotation does not exist")
	if (!is(outdir, "character")) stop("outdir must be a character string")
	if (!is(cores, "numeric")) stop("cores must be a positive integer")
	if (cores < 1 | !cores %% 1 == 0) stop("cores must be a positive integer")
	
	## Make sure the output directory exists.
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Store some of the settings.
	go_obj@settings$genome_assembly <- genome_assembly
	go_obj@settings$genome_annotation <- genome_annotation
	go_obj@settings$star_index <- outdir

	## Print out some information on genome index generation.
	message(
		"\n## Generating STAR Genome Index\n",
		"##\n",
		"## Assembly: ", genome_assembly, "\n",
		"## Annotation: ", genome_annotation, "\n",
		"## Output Directory: ", outdir, "\n",
		"## Cores: ", cores, "\n\n",
		"...Started indexing genome\n"
	)

	## Generate the STAR genome index.
	command <- paste(
		"STAR", "--runMode genomeGenerate",
		"--runThreadN", cores,
		"--genomeDir", outdir,
		"--genomeFastaFiles", genome_assembly,
		"--sjdbGTFfile", genome_annotation
	)
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

	message("...Done indexing genome!")
	return(go_obj)
}

#' Align Reads
#'
#' @description
#' Align reads to genome using STAR
#'
#' @import tibble
#' @importFrom purrr pwalk
#'
#' @param go_obj gostripes object
#' @param outdir Directory to save aligned files to
#' @param cores Number of CPU cores/threads to use
#'
#' @return gostripes object and coordinate sorted and indexed BAMs
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
#'   align_reads("./scratch/aligned")
#'
#' @rdname align_reads-function
#'
#' @export

align_reads <- function(go_obj, outdir, cores = 1) {

	## Ensure validity of input.
	if (!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if (!is(outdir, "character")) stop("outdir must be a character string")
	if (!is(cores, "numeric")) stop("cores must be a positive integer")
	if (cores < 1 | !cores %% 1 == 0) stop("cores must be a positive integer")

	## Make sure the output directory exists.
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Store some of the settings.
	go_obj@settings$star_aligned <- outdir

	## Print out some information on read alignment.
	message(
		"\n## Read Alignment\n",
		"##\n",
		"## Output Directory: ", outdir, "\n",
		"## Cores: ", cores, "\n"
	)

	## Align reads using STAR.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)
		message("...Processing ", args$sample)

		# Get sequencing mode for sample.
		seq_mode <- args$seq_mode

		# Get names of cleaned and dusted fastq files.
		if (seq_mode == "paired") {
			cleaned_R1 <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, "_READ1.fq"))
			cleaned_R2 <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, "_READ2.fq"))
			input_reads <- paste(cleaned_R1, cleaned_R2)
		} else {
			input_reads <- file.path(go_obj@settings$fastq_outdir, paste0("decon_", args$sample_name, ".fq"))
		}

		# Align reads using STAR.
		command <- paste(
			"STAR",
			"--runThreadN", cores,
			"--genomeDir", go_obj@settings$star_index,
			"--readFilesIn", input_reads,
			"--outSAMtype BAM SortedByCoordinate",
			"--outFileNamePrefix", file.path(outdir, paste0(args$sample_name, "_"))
		)
		message("......Aligning reads")
		system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)

		# Index the aligned bams.
		message("......Indexing coordinate sorted BAMs")
		command <- paste(
			"samtools index",
			file.path(outdir, paste0(args$sample_name, "_Aligned.sortedByCoord.out.bam"))
		)
		system(command)
		message("......Done aligning and indexing ", args$sample_name, "!")
	})

	return(go_obj)
}
