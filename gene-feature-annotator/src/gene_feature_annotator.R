#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(ChIPseeker)
  library(rtracklayer)
})

# Get input arguments
args <- commandArgs(trailingOnly=TRUE)

# Check if correct number of arguments provided
if (length(args) != 3) {
  stop("Please provide paths to the .bed, .gtf file and the output .tsv files.")
}

# Get file paths
bed_file <- args[1]
gtf_file <- args[2]
output_file <- args[3]

print("Create a TxDb from .gtf file")
gtf_features <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format='gtf')

print("Load .bed file")
bed_granges <- ChIPseeker::readPeakFile(bed_file)

print("Annotate peaks with gene features")
annotated_peaks <- ChIPseeker::annotatePeak(
  bed_granges,
  TxDb=gtf_features,
  tssRegion=c(-3000, 3000),
  level='transcript',
  verbose=FALSE
)

print("Extract gene names for each gene id")
gtf.df <- as.data.frame(rtracklayer::import(gtf_file))
genes <- unique(gtf.df[ ,c("gene_id", "gene_name", "gene_type")])

print("Merge resulting table with gene names")
annotated_peaks_df <- data.frame(annotated_peaks@anno)
annotated_peaks_df <- merge(
  annotated_peaks_df,
  genes,
  by.x="geneId",
  by.y="gene_id",
  all.x=TRUE,
  sort=FALSE
)

print("Save annotation table to output")
write.table(
  annotated_peaks_df,
  file=output_file,
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

print(annotated_peaks)
