# setwd to oligos_redux

library(Biostrings)
library(dplyr)
library(stringr)

# read the mapping data
mapping_data <- read.csv("consensus/highest_median_kmer_samples.csv")

# pull out fastas only
parent_dir <- "results"
all_files <- list.files(parent_dir, pattern = "\\.(fasta|fa)$", full.names = TRUE, recursive = TRUE)

# Initialize a list to store processed sequences
combined_sequences <- list() # Using a standard list for flexibility before unlisting

# pull out name and kmer sizes
for (i in 1:nrow(mapping_data)) {
  sample <- mapping_data$SampleType[i]
  kmer <- mapping_data$KmerSize[i]
  
  # in the directory find the file that matches the mapping file for highest median kmer samples
  match_pattern <- paste0("(?=.*", sample, ")(?=.*", kmer, ")")
  target_file <- all_files[str_detect(all_files, match_pattern)]
  
  if (length(target_file) > 0) {
    # Process the first matching file
    f <- target_file[1]
    dna <- readDNAStringSet(f)
    
    
    # pull the rank out of the headers, to be more compatible with blastn
    original_headers <- names(dna)
    
    # Extract digits following "rank="
    rank_val <- str_extract(original_headers, "(?<=rank=)\\d+")
    
    # create unique names: [Sample]_[Kmer]_[Rank]_[IterativeIndex]
    seq_idx <- sprintf("%02d", 1:length(dna))
    names(dna) <- paste0(sample, "_k", kmer, "_r", rank_val, "_", seq_idx)
    
    # --------------------------

    # filter out dna where the alphabet is noncanonical to DNA (good for Nones)
    n_counts <- alphabetFrequency(dna)[, "N"]
    non_empty_indices <- which(n_counts < width(dna))
    dna_filtered <- dna[non_empty_indices]

    # 5. Add to master list if valid sequences remain
    if (length(dna_filtered) > 0) {
      combined_sequences[[i]] <- dna_filtered
      message(paste("Processed:", basename(f), "| Extracted Rank:", rank_val))
    }
  } else {
    warning(paste("No match found for Sample:", sample, "Kmer:", kmer))
  }
}

# 6. Combine and save
combined_sequences <- combined_sequences[!sapply(combined_sequences, is.null)]
combined_longest_consensus <- unlist(DNAStringSetList(combined_sequences))

storage_dir <- "consensus"
if (!dir.exists(storage_dir)) dir.create(storage_dir, recursive = TRUE)

output_path <- file.path(storage_dir, "combined_consensus_sequences.fasta")

# make the fasta sequenes all on one line by making the width huge
writeXStringSet(combined_longest_consensus, output_path, width = 20000)

message(paste("Final combined FASTA saved successfully to:", output_path))