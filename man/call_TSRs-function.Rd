% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tss_analysis.R
\name{call_TSRs}
\alias{call_TSRs}
\title{Call TSRs}
\usage{
call_TSRs(go_obj, threshold, clust_dist)
}
\arguments{
\item{go_obj}{gostripes object}

\item{threshold}{TSSs with read count below threshold will be discarded}

\item{clust_dist}{TSSs within this number of base pairs will be clustered}
}
\value{
gostripes object with list of TSR GRanges
}
\description{
Basic TSS clustering into TSRs based on naive global read threshold.
}
\details{
This function will first remove TSSs with an aggregate read number below 'threshold'.
Surviving TSSs will then be clustered if they are within 'clust_dist' bases of eachother.
}
\examples{
R1_fastq <- system.file("extdata", "S288C_R1.fastq", package = "gostripes")
R2_fastq <- system.file("extdata", "S288C_R2.fastq", package = "gostripes")
rRNA <- system.file("extdata", "Sc_rRNA.fasta", package = "gostripes")
assembly <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa", package = "gostripes")
annotation <- system.file("extdata", "Saccharomyces_cerevisiae.R64-1-1.99.gtf", package = "gostripes")

sample_sheet <- tibble::tibble(
  "sample_name" = "stripeseq", "replicate_ID" = 1,
  "R1_read" = R1_fastq, "R2_read" = R2_fastq
)

go_object <- gostripes(sample_sheet) \%>\%
  process_reads("./scratch/cleaned_fastq", rRNA) \%>\%
  fastq_quality("./scratch/fastqc_reports") \%>\%
  genome_index(assembly, annotation, "./scratch/genome_index") \%>\%
  align_reads("./scratch/aligned") \%>\%
  process_bams("./scratch/cleaned_bams") \%>\%
  call_TSSs \%>\%
  call_TSRs(threshold = 3, clust_dist = 25)

}
