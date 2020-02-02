# gostripesR v0.2.0

Processing and quality control of STRIPE-seq FASTQ files.

## Preparing Required Software

### Singularity Container

It is recommended to use the provided singularity container with all required software installed.
Singularity are containers similar to docker containers that allow compatability and reproducibilty for software and workflows.
You must first install the [singularity software](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps) 
onto your machine to use containers.

When you have singularity isntalled and are ready to run the workflow,
you can then download the gostripes container to access all required software.
Create a new scratch directory and navigate into it, and then follow the instructions below.

Pull the singularity container from Sylabs Cloud.
```
singularity pull --arch amd64 library://rpolicastro/default/gostripes:0.2.0
```

Start R within the container to gain access to the installed software.
```
singularity exec -eCB "$(pwd)" -H "$(pwd)" gostripes_0.2.0.sif
```

You are now ready to use gostripes!

## Quickstart

```
library("gostripes")

## Load example data from package.

R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")

rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")

assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")


## Create example sample sheet.

sample_sheet <- tibble::tibble(
	"sample_name" = "stripeseq",
	"replicate_ID" = 1,
	"R1_read" = R1_fastq,
	"R2_read" = R2_fastq
)

## Running the workflow on the example data.

go_object <- gostripes(sample_sheet) %>%
	process_reads("./scratch/cleaned_fastq", rRNA, cores = 4) %>%
	fastq_quality("./scratch/fastqc_reports", cores = 4) %>%
	genome_index(assembly, annotation, "./scratch/genome_index", cores = 4) %>%
	align_reads("./scratch/aligned", cores = 4) %>%
	process_bams("./scratch/cleaned_bams", cores = 4) %>%
	count_features(annotation, cores = 4) %>%
	export_counts("./scratch/counts") %>%
	call_TSSs %>%
	export_TSSs("./scratch/TSSs") %>%
	call_TSRs(3, 25) %>%
	export_TSRs("./scratch/TSRs")
```

## Detailed Start

### Preparing Data

gostripes takes demultiplexed STRIPE-seq FASTQ files as input, in either paired or single-end sequencing format.
For paired end data it is important that the forward and reverse reads are in the same order in both files.

gostripes is also able to handle multiple samples at the same time using a sample sheet.
The sample sheet should have 4 columns: sample_name, replicate_ID, R1_read, R2_read.
Each sample in a group of biological replicates should have the same replicate ID.
The R1 and R2 reads should include the path to the file as well as the file name.
If the samples were sequenced in single-end mode, you can leave the entries in R2_read blank.

```
library("gostripes")

R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")

sample_sheet <- tibble::tibble(
        "sample_name" = "stripeseq",
        "replicate_ID" = 1,
        "R1_read" = R1_fastq,
        "R2_read" = R2_fastq
)

go_object <- gostripes(sample_sheet)
```
### Quality Control of FASTQ Files

The first main step of STRIPE-seq analysis is the quality control and filtering of the FASTQ files.
First, R1 read structure is ensured by looking for 'NNNNNNNNTATAGGG' at the beginning of the R1 read,
which corresponds to the UMI:spacer:riboG of the template switching oligonucleotide.
Second, the UMI is stashed in the read name, allowing it to be used for duplicate removal in single-end data (and optionally paired-end).
Third, the remaining TATAGGG after UMI removal is trimmed.
Finally, contaminant reads such as rRNA are filtered out.
This requires a FASTA file containing the contaminant sequences to search against.

As further quality assurance FastQC quality reports are generated both for the raw FASTQ files,
and the processed FASTQ files.

```
rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")

go_object <- process_reads(go_object, "./scratch/cleaned_fastq", rRNA, cores = 4)
go_object <- fastq_quality(go_object, "./scratch/fastqc_reports", cores = 4)
```

### Aligning Reads to Genome

After quality control of the FASTQ files, the reads can then be mapped to the genome.
First, a STAR genome index is generated from the FASTA genome assembly and GTF genome annotation file.
Then, the FASTQ files are mapped to the genome using this index.

```
assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")

go_object <- genome_index(go_object, assembly, annotation, "./scratch/genome_index", cores = 4)
go_object <- align_reads(go_object, "./scratch/aligned", cores = 4)
```

### Quality Control of BAM Files

After aligning the reads to the genome, the result is a coordinate sorted and indexed BAM.
Two main quality control steps are taken with this BAM to ensure the most accurate measurment of true TSSs.
First, PCR duplicates reads are removed either using samtools (paired-end) or UMI-tools (single-end).
During this step, various other checks are made, such as ensuring properly paired reads and removing non-primary alignments.
Second, any TSS that has more than 3 soft-clipped bases adjacent to it is removed from the BAM.

```
go_object <- process_bams(go_object, "./scratch/cleaned_bams", cores = 4)
```
### Feature Counting

After the quality contol steps, the resulting BAMs can be used for RNA-seq like feature counting.
Each read or read-pair will be assigned to the closest overlapping exon,
and a summary of overlapping read counts will be produced for each gene.
These feature counts can then optionally be exported as a table.

```
go_object <- count_features(go_object, annotation, cores = 4)
export_counts(go_object, "./scratch/counts")
```

### Rudimantary TSS and TSR Calling

The final BAMs are also ready for TSS and TSS cluster (TSR or cTSS) analysis.
There are many great software suites available for this, including
[TSRchitect](https://bioconductor.org/packages/release/bioc/html/TSRchitect.html),
[CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html),
[ADAPT-CAGE](https://gitlab.com/dianalab/adapt-cage), and
[CAGEfightR](https://bioconductor.org/packages/release/bioc/html/CAGEfightR.html).
For convenience, gostripes includes some rudimentary functions for basic TSS and TSR calling.

Although 5' ends with 3 or less soft-clipped bases are retained in the bam quality control steps, those bases are not considered when calling TSSs.
For TSR calling, TSSs with less than the user defined threshold number of reads are first removed.
Surviving TSSs within the user defined number of bases (25 by default) are then clustered into a TSRs/cTSS.
The resulting TSSs and TSRs/cTSSs can be exported as BEDGRAPH and BED files respectively.

```
go_object <- call_TSSs(go_object)
export_TSSs(go_object, "./scratch/TSSs")

go_object <- call_TSRs(go_object, 3, 25)
export_TSRs(go_object, "./scratch/TSRs")
```

## Acknowledgments

The development of gostripes would not be possible without these great software packages.

* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/): FASTQ qauality control.
* [TagDust 2](http://tagdust.sourceforge.net/): FASTQ read filtering.
* [STAR](https://github.com/alexdobin/STAR): Short read sequence aligner.
* [Samtools](http://www.htslib.org/): SAM/BAM file manipulation.
* [Picard](https://broadinstitute.github.io/picard/): Manipulation of SAM/BAM files.

A special shoutout to the [tidyverse](https://www.tidyverse.org/) for making data science in R easy.
Also, a sincere thank you to [Bioconductor](http://bioconductor.org/) and it's varied contributors for hosting so many invaluable tools.
