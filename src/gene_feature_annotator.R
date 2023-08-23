#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(ChIPseeker)
  library(rtracklayer)
  library(optparse)
})


option_list <- list(
  make_option(
    c("-g", "--gtf"),
    help="Path to the gtf file."
  ),
  make_option(
    c("-i","--filelist"),
    help="List of comma-separated input BED file paths to annotate."
  ),
  make_option(
    c("-o","--dirOutput"),
    help="Path to the output directory."
  )
)

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

gtfpath <- opt$gtf
myfilelist <- strsplit(opt$filelist, ",")[[1]]
dir_out <- opt$dirOutput

print("Create a TxDb from .gtf file")
gtf_features <- GenomicFeatures::makeTxDbFromGFF(gtfpath, format='gtf')

print("Extract gene names for each gene id")
gtf.df <- as.data.frame(rtracklayer::import(gtfpath))
genes <- unique(gtf.df[,c("gene_id", "gene_name", "gene_type")])

for (bed_file in myfilelist) {
  print(bed_file)

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

  output_file <- paste0(dir_out, sub(".bed", "", basename(bed_file)), ".anno.tsv")
  print(output_file)

  print("Save annotation table to output")
  write.table(
    annotated_peaks_df,
    file=output_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

  print(annotated_peaks)
}