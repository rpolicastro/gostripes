# gostripes

Processing of raw STRIPE-seq files.

## Quickstart

```
library("gostripes")

R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")

sample_sheet <- tibble(
	"sample_name" = "stripeseq",
	"replicate_ID" = 1,
	"R1_read" = R1_fastq,
	"R2_read" = R2_fastq
)

go_object <- gostripes(sample_sheet)
```
