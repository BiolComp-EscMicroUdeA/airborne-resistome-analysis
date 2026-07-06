
# Script: Data loading and preprocessing
# Description: Loads RGI bwt gene_mapping_data output files for
#              all samples, applies quality filters, normalizes
#              read counts (TPM/RPKM), and integrates sample metadata.
# Input:  - RGI bwt output files (*_gene_mapping_data.txt) and (*_overall_mapping_stats.txt) per sample
#         - Sample metadata file (Metadata_PM.csv)
# Output: - all_data: normalized AMR gene data for all samples
#         - all_data_metadata: all_data joined with sample metadata


# Load libraries
library(tidyverse)
library(fs)
library(scales)
library(pheatmap)
library(ggpubr)
library(vegan)
library(circlize)
library(Hmisc)
library(rstatix)
library(igraph)
library(UpSetR)
library(VennDiagram)


# --- Paths (adjust according to your environment) ---
base_dir      <- "path/to/rgi/output/folders"   # directory containing one folder per sample
metadata_path <- "path/to/metadata/muestras.csv" # sample metadata file


# --- Helper function: extract total mapped reads from RGI stats file ---
extract_mapped_reads <- function(stats_path) {
  if (!file_exists(stats_path)) return(NA_real_)
  
  stats_content <- readLines(stats_path, warn = FALSE)
  mapped_line   <- grep("^Mapped reads", stats_content, value = TRUE)
  
  if (length(mapped_line) == 0) return(NA_real_)
  
  as.numeric(str_extract(mapped_line, "\\d+"))
}


# --- Main function: process a single sample directory ---
# Reads gene_mapping_data, applies filters, and computes TPM and RPKM
process_sample <- function(sample_dir) {
  sample_name <- path_file(sample_dir)
  
  # Locate required files
  gene_file  <- dir_ls(sample_dir, regexp = "gene_mapping_data\\.txt$")
  stats_file <- dir_ls(sample_dir, regexp = "overall_mapping_stats\\.txt$")
  
  # Validate that exactly one of each file exists
  if (length(gene_file) != 1 || length(stats_file) != 1) {
    message(sprintf("Sample %s: missing or duplicate files — skipping", sample_name))
    return(NULL)
  }
  
  # Extract total mapped reads for normalization
  total_mapped <- extract_mapped_reads(stats_file)
  if (is.na(total_mapped)) {
    message(sprintf("Sample %s: 'Mapped reads' not found in stats file — skipping", sample_name))
    return(NULL)
  }
  
  # Read gene-level RGI output
  gene_arg <- read_tsv(gene_file, show_col_types = FALSE)
  
  # --- Filters ---
  # Keep only: CARD database entries, protein homolog models,
  # and genes with average percent coverage > 70%
  gene_arg <- gene_arg %>%
    filter(
      `Reference DB`             == "CARD",
      `Reference Model Type`     == "protein homolog model",
      `Average Percent Coverage`  > 70
    )
  
  # Ensure Reference Length is numeric
  gene_arg$`Reference Length` <- as.numeric(gene_arg$`Reference Length`)
  
  # --- Normalization: TPM and RPKM ---
  gene_arg <- gene_arg %>%
    mutate(
      Sample = sample_name,
      # Reads Per Kilobase (RPK) — intermediate value for TPM
      RPK  = `All Mapped Reads` / (`Reference Length` / 1000),
      # Transcripts Per Million (TPM)
      TPM  = (RPK / sum(RPK)) * 1e6,
      # Reads Per Kilobase per Million mapped reads (RPKM)
      RPKM = (`All Mapped Reads` * 1e9) / (`Reference Length` * total_mapped),
      .before = 1
    ) %>%
    select(
      Sample, `ARO Term`, `ARO Accession`, `Reference Model Type`,
      `Average Percent Coverage`, `Resistomes & Variants: Observed Pathogen(s)`,
      `Drug Class`, `AMR Gene Family`, `Reference Length`,
      `All Mapped Reads`, RPKM, TPM, `Resistance Mechanism`
    )
  
  return(gene_arg)
}


# --- Process all samples ---
sample_dirs <- dir_ls(base_dir, type = "directory")

# Combine all sample results into a single data frame
all_data <- map_df(sample_dirs, process_sample, .progress = TRUE)


# --- Load and integrate sample metadata ---
metadatos_pm <- read.csv2(metadata_path)

all_data_metadata <- left_join(all_data, metadatos_pm, by = "Sample")