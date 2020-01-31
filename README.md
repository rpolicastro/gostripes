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
	process_reads("./scratch/cleaned_fastq", rRNA) %>%
	genome_index(assembly, annotation, "./scratch/genome_index", cores = 4) %>%
	align_reads("./scratch/aligned", cores = 4) %>%
	process_bams("./scratch/cleaned_bams", cores = 4)

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
This requires a fasta file containing the contaminant sequences to search against.

```
rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
go_object <- process_reads(go_object, "./scratch/cleaned_fastq", rRNA)
```

### Aligning Reads to Genome

After quality control of the fastq files, the reads can then be mapped to the genome.
First, a STAR genome index is generated from the fasta genome alignment and genome annotation file.
Then, the fastq files are mapped to the genome using this index.

```
assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")

go_object <- genome_index(go_object, assembly, annotation, "./scratch/genome_index", cores = 4)
go_object <- align_reads(go_object, "./scratch/aligned", cores = 4)
```

### Quality Control of BAM Files

After aligning the reads to the genome, the result is a coordinate sorted bam.
Two main quality control steps are taken with this bam to ensure the most accurate measurment of true TSSs.
First, PCR duplicates reads are removed either using samtools (paired-end) or UMI-tools (single-end).
Second, any TSS that has more than 3 soft-clipped bases adjacent to it is removed from the bam.

```
go_object <- process_bams(go_object, "./scratch/cleaned_bams", cores = 4)
```

### Rudimantary TSS and TSR calling.

After the quality contol steps, the resulting bams are ready for TSS and TSS cluster (TSR or CTSSs) analysis.
There are many great software suites available for this, including TSRchitect, CAGEr, ADAPT-CAGE, and CAGEfightR.
For convenience, gostripes includes some rudimentary functions for basic TSS and TSR calling from the bam also.

Although 5' ends with 3 or less soft-clipped bases are retained in the bam quality control steps, those bases are not considered when calling TSSs.
For TSR calling, TSSs with less than the user defined threshold number of reads are first removed.
Surviving TSSs within the user defined number of bases (25 by default) are then clustered into a TSR/CTSS.
The resulting TSSs and TSRs/CTSSs can be exported as bedgraphs and bed files respectively.

```
go_object <- call_TSSs(go_object)
export_TSSs(go_object, "./scratch/TSSs")

go_object <- call_TSRs(go_object, 3, 25)
export_TSRs(go_object, "./scratch/TSRs")
```
