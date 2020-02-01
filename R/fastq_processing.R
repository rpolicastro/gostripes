
#' Process Reads
#'
#' Remove rRNA contamination, low complexity reads, and trim UMI-Spacer-GGG
#'
#' @import tibble
#' @importFrom stringr str_replace
#' @importFrom purrr pwalk
#' @importFrom magrittr %>%
#'
#' @param go_obj gostripes object
#' @param contamination_fasta fasta file containing contaminants to remove, such as rRNA
#' @param outdir output directory for filtered and trimmed reads
#' @param cores Number of CPU core/threads to use
#'
#' @rdname process_reads-function
#'
#' @export

process_reads <- function(go_obj, outdir, contamination_fasta, cores = 1) {

	## Check input validity.
	if(!is(go_obj, "gostripes")) stop("go_obj must be a gostripes object")
	if(!is(outdir, "character")) stop("output directory must be a character string")
	if(!is(contamination_fasta, "character")) stop("contamination_fasta must be a character string")
	if(!file.exists(contamination_fasta)) stop("contamination_fasta file not found")
	if(!is(cores, "numeric")) stop("cores must be a positive integer")
	if(!cores %% 1 == 0) stop("cores must be a positive integer")
	
	## Add output directory to settings and make sure it exists.
	go_obj@settings$fastq_outdir <- outdir
	if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

	## Print out some analysis information.
	message(
		"\n## FASTQ Processing\n",
		"##\n",
		"## Output Directory: ", outdir, "\n",
		"## Contaminant FASTA: ", contamination_fasta, "\n",
		"## Cores: ", cores, "\n",
		"##\n",
		"## Starting Analysis...\n"
	)

	## Process each pair of reads.
	pwalk(go_obj@sample_sheet, function(...) {
		args <- list(...)
		message("...Analyzing ", args$sample_name)
		
		# Check for paired end status.
		seq_mode <- args$seq_mode

		# Check for proper R1 read structure.
		message("......Checking for proper R1 read structure")
		R1_read <- args$R1_read
		R2_read <- ifelse(seq_mode == "paired", args$R2_read, NA) 

		read_structure(
			args$R1_read, args$R2_read, 
			args$sample_name, outdir
		)

		# UMI-tools to stash UMI in read name.
		message("......Stashing the UMI sequence in the read name")
		proper_R1 <- file.path(outdir, paste0("proper_", args$sample_name, "_R1.fastq"))
		proper_R2 <- ifelse(
			seq_mode == "paired",
			file.path(outdir, paste0("proper_", args$sample_name, "_R2.fastq")), NA
		)

		stash_umi(
			proper_R1, proper_R2,
			args$sample_name, outdir
		)

		# Remove spacer and ribo-Gs.
		message("......Trimming the spacer (TATA) and ribo-Gs (GGG)")
		stashed_R1 <- file.path(outdir, paste0("stashed_", args$sample_name, "_R1.fastq"))
		stashed_R2 <- ifelse(
			seq_mode == "paired",
			file.path(outdir, paste0("stashed_", args$sample_name, "_R2.fastq")), NA
		)

		remove_extra(
			stashed_R1, stashed_R2,
			args$sample_name, outdir
		)

		# TagDust2 to remove contaminants such as rRNA and low complexity reads.
		message("......Removing contaminant reads")
		trimmed_R1 <- file.path(outdir, paste0("trimmed_", args$sample_name, "_R1.fastq"))
		trimmed_R2 <- ifelse(
			seq_mode == "paired",
			file.path(outdir, paste0("trimmed_", args$sample_name, "_R2.fastq")), NA
		)

		remove_contaminants(
			trimmed_R1, trimmed_R2, args$sample_name,
			contamination_fasta, outdir, cores
		)

		message("......Done processing ", args$sample_name, "!")

	})

	return(go_obj)
}

#' Check R1 Structure
#'
#' Ensure proper structure of R1 read, and remove spacer and ribo-Gs
#'
#' @import tibble
#' @importFrom stringr str_which
#' @importFrom Biostrings readDNAStringSet DNAStringSet writeXStringSet
#'
#' @param R1_read File name and directory for R1 read
#' @param R2_read File name and directory for R2 read
#' @param sample_name Name of the sample
#' @param outdir Output directory for the properly structured reads
#' @param structure_regex Regular expression for structure of read 5' end
#'
#' @rdname read_structure-function
#'
#' @export

read_structure <- function(R1_read, R2_read, sample_name, outdir, structure_regex = "^[ATGCN]{8}TATAGGG") {

	## Load up R1 read and check for reads matching proper R1 structure.
	R1_data <- readDNAStringSet(R1_read, format = "fastq", with.qualities = TRUE)
	R1_keep_index <- str_which(R1_data, structure_regex)
	R1_keep <- R1_data[R1_keep_index]

	## If paired end keep only R2 reads where the R1 pair wasn't filtered.
	if (!is.na(R2_read)) {
		R2_data <- readDNAStringSet(R2_read, format = "fastq", with.qualities = TRUE)
		R2_keep <- R2_data[R1_keep_index]
	}

	## New sample names.
	proper_R1 <- file.path(outdir, paste0("proper_", sample_name, "_R1.fastq"))
	if (!is.na(R2_read)) {
		proper_R2 <- file.path(outdir, paste0("proper_", sample_name, "_R2.fastq"))
	}

	## Write samples to fastq files.
	writeXStringSet(R1_keep, proper_R1, format = "fastq")
	if (!is.na(R2_read)) {
		writeXStringSet(R2_keep, proper_R2, format = "fastq")
	}
}

#' Remove Contaminants
#'
#' Remove low complexity reads and contaminant reads such as rRNA
#'
#' @importFrom stringr str_replace
#'
#' @param proper_R1 R1 reads with proper structure
#' @param proper_R2 R2 reads where R1 had proper structure
#' @param sample_name Name of the sample
#' @param contamination_fasta Fasta file containing the contaminants
#' @param outdir Output directory for decontaminated files
#' @param cores Number of CPU cores available
#'
#' @rdname remove_contaminants-function
#'
#' @export

remove_contaminants <- function(trimmed_R1, trimmed_R2, sample_name, contamination_fasta, outdir, cores) {

	## Build TagDust2 command.
	command <- paste(
		"tagdust",
		"-ref", contamination_fasta,
		"-fe 3", "-t", cores, "-dust 97",
		"-o", file.path(outdir, paste0("decon_", sample_name)), "-1 R:N",
		trimmed_R1
	)

	## Add paired end read if available.
	if (!is.na(trimmed_R2)) {command <- paste(command, trimmed_R2)}
	
	## Run TagDust2 command.
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

#' Stash UMI
#'
#' Add UMI to read name to facilitate PCR duplicate removal
#'
#' @param deconned_R1 R1 read with contaminants removed
#' @param deconned_R2 R2 read with contaminants removed
#' @param sample_name Name of the sample
#' @param outdir Output directory
#' @param umi_pattern UMI pattern on 5' end of read
#'
#' @rdname stash_umi-function
#'
#' @export

stash_umi <- function(proper_R1, proper_R2, sample_name, outdir, umi_pattern = "NNNNNNNN") {

	## Output files.
	stashed_R1 <- file.path(outdir, paste0("stashed_", sample_name, "_R1.fastq"))
	if (!is.na(proper_R2)) {
		stashed_R2 <- file.path(outdir, paste0("stashed_", sample_name, "_R2.fastq"))
	}

	## Build UMI-tools command.
	command <- paste(
                        "umi_tools extract",
                        "--extract-method=string",
                        paste0("--bc-pattern=", umi_pattern),
			"-L", file.path(outdir, paste0(sample_name, "_umi.log")),
			"-I", proper_R1, "-S", stashed_R1
	)

	## Add paired end info to command if available.
	if (!is.na(proper_R2)) {
		command <- paste(
			command,
			paste0("--read2-in=", proper_R2),
			paste0("--read2-out=", stashed_R2)
		)
	}

	## Run UMI-tools command.
	system(command, ignore.stdout = TRUE, ignore.stderr = TRUE)	
}

#' Remove Extra
#'
#' Remove spacer (TATA) and ribo-Gs (GGG) from 5' end of read
#'
#' @importFrom ShortRead readFastq writeFastq
#' @importFrom IRanges narrow
#' @importFrom magrittr %>%
#'
#' @param stashed_R1 R1 read with UMI stashed in read name
#' @param stashed_R2 R2 read with UMI stashed in read name
#' @param sample_name Name of the sample
#' @param outdir Output directory of reads with spacer adn ribo-Gs trimmed
#' @param extra_bases The structure of the extra 5' bases
#'
#' @rdname remove_extra-function
#'
#' @export

remove_extra <- function(stashed_R1, stashed_R2, sample_name, outdir, extra_bases = "TATAGGG") {

	## Get trim ammount.
	trim_length <- nchar(extra_bases) + 1

	## Make new file names.
	trimmed_R1 <- file.path(outdir, paste0("trimmed_", sample_name, "_R1.fastq"))
	if (!is.na(stashed_R2)) {
		trimmed_R2 <- file.path(outdir, paste0("trimmed_", sample_name, "_R2.fastq"))
	}

	## Trim the R1 read.
	stashed_R1 %>%
		readFastq %>%
		narrow(start = trim_length) %>%
		writeFastq(trimmed_R1, compress = FALSE)

	## Read in R2 read if available.
	if (!is.na(stashed_R2)) {
		stashed_R2 %>%
			readFastq %>%
			writeFastq(trimmed_R2, compress = FALSE)
	}
}
