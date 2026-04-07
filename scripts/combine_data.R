#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fs)
})

# ----------------------------
# Command Line Arguments
# ----------------------------
# This grabs every name you type after the script name
experiments <- commandArgs(trailingOnly = TRUE)

if (length(experiments) == 0) {
  stop("Please provide at least one experiment name. 
Usage: Rscript process_oligos.R exp1 exp2 exp3 ...", call. = FALSE)
}

cutoff <- 0.05

# ----------------------------
# Functions
# ----------------------------
combine_data <- function(files, fdr_cutoff, out_prefix = NULL) {
  dfs <- list()
  
  for (k in seq.int(5, 15, 1)) {
    pattern <- paste0('_', k, 'mers_20c_')
    idx <- grep(pattern, files)
    
    if (length(idx) == 0) next
    
    message(paste0('  - Reading k=', k))
    
    df <- read_csv(files[idx], show_col_types = FALSE) %>%
      mutate(kmer_size = k)
    
    bg_cols    <- grep('^bg_\\d*', names(df), value = TRUE)
    zero_cols  <- grep('^0w_\\d*', names(df), value = TRUE)
    one_cols   <- grep('^1w_\\d*', names(df), value = TRUE)
    three_cols <- grep('^3w_\\d*', names(df), value = TRUE)
    
    df <- df %>%
      mutate(
        bg_mean_count = rowMeans(select(., all_of(bg_cols))),
        `0w_mean_count` = rowMeans(select(., all_of(zero_cols))),
        `1w_mean_count` = rowMeans(select(., all_of(one_cols))),
        `3w_mean_count` = rowMeans(select(., all_of(three_cols)))
      ) %>%
      select(-matches('^bg_\\d$'), -matches('^0w_\\d$'),
             -matches('^1w_\\d$'), -matches('^3w_\\d$'), -any_of("confect"))
    
    df <- df %>% filter(fdr <= fdr_cutoff)
    
    if (!is.null(out_prefix)) {
      dir_create(path_dir(out_prefix))
      outname <- paste0(out_prefix, '_', k, 'mers_weights.csv')
      
      df %>%
        select(edge, `3w`) %>%
        rename(weight = `3w`) %>%
        arrange(desc(weight)) %>%
        fwrite(outname)
    }
    
    dfs[[as.character(k)]] <- df
  }
  
  return(bind_rows(dfs))
}

# ----------------------------
# Main Loop
# ----------------------------

# Ensure the combined output folder exists
dir_create('figures/combined')

for (exp_name in experiments) {
  message(paste("\n>>> Processing:", exp_name))
  
  input_path  <- file.path('results', exp_name, 'labeled')
  weight_path <- file.path('results', exp_name, 'consensus_seqs_nolib', exp_name)
  final_csv   <- file.path('figures/combined', paste0(exp_name, '_data.csv'))
  
  # Find files
  files <- dir(input_path, pattern = '*20c_all_washes_enriched_nolib_labeled.csv', full.names = TRUE)
  
  if (length(files) == 0) {
    warning(paste("No matching files found in", input_path, "- Skipping."))
    next
  }
  
  # Process and mutate 'exp' column
  processed_data <- combine_data(files, fdr_cutoff = cutoff, out_prefix = weight_path) %>%
    mutate(exp = exp_name)
  
  # Save
  write.csv(processed_data, final_csv, row.names = FALSE)
  message(paste("Saved to:", final_csv))
}

message("\nAll experiments finished!")