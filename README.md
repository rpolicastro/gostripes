# gostripes

Processing of raw STRIPE-seq files.

## Quickstart

```
library("gostripes")

## Load example fastq files from package.
R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")

## Create example sample sheet.
sample_sheet <- tibble(
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
```
