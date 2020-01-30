# gostripes

Processing of raw STRIPE-seq files.

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

Shell into the container to gain access to the installed software.
```
singularity shell -eCB "$(pwd)" -H "$(pwd)" gostripes_0.2.0.sif
```

Activate the conda environment within the container, and start R.
```
. /opt/conda/etc/profile.d/conda.sh
conda activate gostripes
R
```

You are now ready to use gostripes!

## Quickstart

```
library("gostripes")

## Load example fastq files from package.

R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")

## Create example sample sheet.

sample_sheet <- tibble::tibble(
	"sample_name" = "stripeseq",
	"replicate_ID" = 1,
	"R1_read" = R1_fastq,
	"R2_read" = R2_fastq
)

## Create gostripes object.

go_object <- gostripes(sample_sheet)

## Process fastq files.

rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
go_object <- process_reads(go_object, "./scratch/cleaned_fastq", rRNA)

## Create a STAR genome index.

assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")
go_object <- genome_index(go_object, assembly, annotation, "./scratch/genome_index", cores = 4)

## Align read to genome using STAR.

go_object <- align_reads(go_object, "./scratch/aligned", cores = 4)

## Process bams.

go_object <- process_bams(go_object, "./scratch/cleaned_bams", cores = 4)

## Call TSSs and export as bedgraph.

go_object <- call_TSSs(go_object)
export_TSSs(go_object, "./scratch/TSSs")

## Call TSRs using naive thresholding and export as beds.

go_object <- call_TSRs(go_object, 3, 25)
export_TSRs(go_object, "./scratch/TSRs")
```

## Detailed Start

### Preparing Data

gostripes takes demultiplexed STRIPE-seq fastq files as input, in either paired or single-end sequencing format.
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

The first main step of STRIPE-seq analysis is the quality control and filtering of the fastq files.
First, R1 read structure is ensured by looking for 'NNNNNNNNTATAGGG' at the beginning of the R1 read,
which corresponds to the UMI:spacer:riboG of the template switching oligonucleotide.
Second, the UMI is stashed in the read name, allowing it to be used for duplicate removal in single-end data (and optionally paired-end).
Third, the remaining TATAGGG after UMI removal is trimmed.
Finally, contaminant reads such as rRNA are filtered out.

```
rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
go_object <- process_reads(go_object, "./scratch/cleaned_fastq", rRNA)
```
