#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(parallel)

# check os pbapply is installed
has_pb <- requireNamespace("pbapply", quietly = TRUE)

# Capture arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript kmer_uniques.R <parent_dir> <target_subdir_name>")
}

parent_dir <- args[1]
target_subdir <- args[2]

# Find relevant files and keep track of time because it is a time intensive script
# specifically goes into results/exp/graphs
start_time <- Sys.time()
all_files <- list.files(path = parent_dir, 
                        pattern = "\\.csv$", 
                        recursive = TRUE, 
                        full.names = TRUE)

target_files <- all_files[basename(dirname(all_files)) == target_subdir]

cat("--------------------------------------------------\n")
cat("Search complete in:", round(difftime(Sys.time(), start_time, units="secs"), 2), "s\n")
cat("Total files found to process:", length(target_files), "\n")
cat("--------------------------------------------------\n")

# Identify the files

if (length(target_files) == 0) {
  stop("No files found! Check if '", target_subdir, "' is spelled exactly right.")
}

process_file <- function(f) {
  # Timing each individual file
  f_start <- Sys.time()
  
  # recursively identify which files are needed
  substrate <- basename(dirname(dirname(f)))
  filename <- basename(f) %>% str_remove(".csv$")
  parts <- str_split(filename, "_")[[1]]
  
  k_val <- as.numeric(str_extract(tail(parts, 1), "\\d+(?=mers)")) # pull out the kmer, replicate, wash and combine
  replicate <- parts[length(parts) - 1]
  condition <- parts[length(parts) - 2]
  sample_type <- paste(parts[1:(length(parts) - 3)], collapse = "_")

  # Reading and counting
  dt <- fread(f, select = c("Node1", "Node2"), header = TRUE, sep = ",", showProgress = FALSE) #rebuild the kmers from the debruijn graphs using the different nodes
  unique_count <- dt[, .(kmer = paste0(Node1, substr(Node2, nchar(Node2), nchar(Node2))))][, uniqueN(kmer)] # count up all of the unique kmers made
  theo_max <- 4^k_val # max possible value of unique kmers based on value of k

  f_end <- Sys.time()
  duration <- round(difftime(f_end, f_start, units="secs"), 2)
  
  # Return result
  return(list(Substrate=substrate, Sample=sample_type, Cond=condition, 
              Rep=replicate, K=k_val, Unique=unique_count, 
              Theo=theo_max, Prop=unique_count/theo_max,
              Proc_Time_Sec = duration))
}

# using the number of cores for the job, use one core as default if not specified
num_cores <- as.numeric(Sys.getenv("NUM_JOBS", unset = 1))
cat("Processing using", num_cores, "cores...\n")

# making progress bar and parallelizing
if (has_pb && .Platform$OS.type != "windows") {
  # Progress bar version (Linux/Mac)
  results_list <- pbapply::pblapply(target_files, process_file, cl = num_cores)
} else {
  # Standard parallel version
  results_list <- mclapply(target_files, process_file, mc.cores = num_cores)
}

# 3. Combine and Save
results <- rbindlist(results_list)

# Final formatting and output
dir.create("data/kmer_space", recursive = TRUE, showWarnings = FALSE)
fwrite(results, "data/kmer_space/kmer_unique_report.csv")

total_duration <- round(difftime(Sys.time(), start_time, units="mins"), 2)
cat("--------------------------------------------------\n")
cat("SUCCESS!\n")
cat("Total runtime:", total_duration, "minutes\n")
cat("Average time per file:", round(mean(results$Proc_Time_Sec), 2), "seconds\n")
cat("Report saved to: data/kmer_space/kmer_unique_report.csv\n")
cat("--------------------------------------------------\n")