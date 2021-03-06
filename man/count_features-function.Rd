% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_count.R
\name{count_features}
\alias{count_features}
\title{Count Features}
\usage{
count_features(go_obj, genome_annotation, cores = 1)
}
\arguments{
\item{go_obj}{gostripes object}

\item{genome_annotation}{Genome annotation in GTF file format}

\item{cores}{Number of CPU cores available}
}
\value{
gostripes object with feature counts matrix
}
\description{
Count reads associated with annotated features
}
\details{
Genome annotations can be found in repositories such as NCBI, UCSC, and ensembl.
The 'genome_annotation' file should be the same GTF used in read alignment for consistency.

This function uses the featureCounts function from Rsubread to summarize counts to annotated features.
First, sequenced fragments are assigned to the nearest exon if there is at least 10 overlapping bases.
If the fragment overlaps more than one feature, it is asigned to the feature with the largest overlap.
Finally, counts for all exons from the same gene are aggregated into a sum score for that gene.
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
  count_features(annotation)

}
