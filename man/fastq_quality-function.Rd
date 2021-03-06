% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastq_quality.R
\name{fastq_quality}
\alias{fastq_quality}
\title{Fastq Quality}
\usage{
fastq_quality(go_obj, outdir, fastq_type = "both", cores = 1)
}
\arguments{
\item{go_obj}{gostripes object}

\item{outdir}{Output directory for FastQC reports}

\item{fastq_type}{Either 'raw', 'processed', or 'both'}

\item{cores}{Number of CPU cores}
}
\value{
FastQC quality reports
}
\description{
Generate quality control reports for FASTQ files using FastQC
}
\details{
The fastq_type specifies the FASTQ files to quality control.
'raw' will analyze the FASTQ files before the initial processing
and quality control steps. 'processed' will look at the FASTQ files
after the quality control and processing steps. 'both' will generate
reports for the 'raw' and 'processed' FASTQ files.
}
\examples{
R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")

sample_sheet <- tibble::tibble(
  "sample_name" = "stripeseq", "replicate_ID" = 1,
  "R1_read" = R1_fastq, "R2_read" = R2_fastq
)

go_object <- gostripes(sample_sheet) \%>\%
  process_reads("./scratch/cleaned_fastq", rRNA) \%>\%
  fastq_quality("./scratch/fastqc_reports")

}
