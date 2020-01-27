
#' Generate Genome Index
#'
#' Generate a STAR genomic index
#'
#' @param go_obj gostripes object
#' @param genome_assembly Genome assembly in fasta file format
#' @param genome_annotation Genome annotation in gtf file format
#' @param outdir Directory to save genome index to
#' @param cores Number of CPU cores/threads to use
#'
#' @rdname genome_index-function
#'
#' @export

genome_index <- function(go_obj, genome_assembly, genome_annotation, outdir, cores = 1) {
	
	## Store some of the settings.
	go_obj@settings$genome_assembly <- genome_assembly
	go_obj@settings$genome_annotation <- genome_annotation
	go_obj@settings$star_index <- outdir

	## Generate the STAR genome index.
	command <- paste(
		"STAR", "--runMode genomeGenerate",
		"--runThreadN", cores,
		"--genomeDir", outdir,
		"--genomeFastaFiles", genome_assembly,
		"--sjdbGTFfile", genome_annotation
	)
	system(command)
}
