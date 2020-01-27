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
```
