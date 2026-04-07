#!/usr/bin/env Rscript

library(Biostrings)
library(dplyr)
library(stringr)
library(readr)
library(fs)

# get the experiment names
experiments <- commandArgs(trailingOnly = TRUE)

if (length(experiments) == 0) {
  stop("Please provide at least one experiment name.", call. = FALSE)
}

# point to results for median kmer counting
parent_folder <- "results"

# find the fasta files in that result folder
all_files <- list.files(path = parent_folder, 
                        pattern = "\\.(fasta|fa)$", 
                        full.names = TRUE, 
                        recursive = TRUE)

# pull out fastas only for the specified substrates
file_list <- all_files[grep(paste(experiments, collapse="|"), all_files)]

# report error if no fasta file
if (length(file_list) == 0) {
  stop("No fasta files found for the provided experiments.")
}

message(paste("Processing experiments:", paste(experiments, collapse = ", ")))

# combine the data to make it easier to read and pull out the DNA sequence from the fasta
all_data <- do.call(rbind, lapply(file_list, function(f) {
  dna <- readDNAStringSet(f)
  
  # read the file name experiment and kmer size
  sample_name <- experiments[sapply(experiments, grepl, x = f)][1]
  kmer_size <- as.numeric(str_extract(basename(f), "\\d+(?=mers)"))
  
  data.frame(
    SampleType = sample_name,
    KmerSize = kmer_size,
    SeqLength = width(dna)
  )
}))

# identify the max median length in the experiment, kmer combinations
summary_report <- all_data %>%
  group_by(SampleType, KmerSize) %>%
  summarize(MedianLength = median(SeqLength, na.rm = TRUE), .groups = "drop") %>%
  group_by(SampleType) %>%
  slice_max(order_by = MedianLength, n = 1, with_ties = FALSE) %>%
  ungroup()

# 7. Save results
output_dir <- "consensus"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "highest_median_kmer_samples.csv")

# We use write_csv without append so it starts fresh each time the shell script runs
write_csv(summary_report, output_file)

message(paste("Summary report written to:", output_file))
print(summary_report)