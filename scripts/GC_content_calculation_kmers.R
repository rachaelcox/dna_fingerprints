# 1. Install/Load Required Packages
if (!require("data.table")) install.packages("data.table", repos = "https://cloud.r-project.org")
if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!require("stringr")) install.packages("stringr", repos = "https://cloud.r-project.org")
if (!require("grid")) install.packages("grid", repos = "https://cloud.r-project.org")

# read in the libraries
library(data.table)
library(ggplot2)
library(stringr)
library(grid)

# provide file paths
results_path  <- "results"
output_folder <- "consensus/GC_content"

# does the directory exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# locate input files
files <- list.files(
  path = results_path,
  pattern = "consensus_seqs\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# confirm the reading of the correcr files
message("Number of files found: ", length(files))

if (length(files) == 0) {
  stop("No files matched pattern 'consensus_seqs\\.csv$' in: ", results_path)
}

# 4. Process files
all_data_list <- list()

# extract the k-mer values from the file name
for (f in files) {
  fname  <- basename(f)
  s_name <- str_split(fname, "_")[[1]][1]
  k_val  <- as.numeric(str_extract(fname, "\\d+(?=mer)"))

  message("Processing: ", fname)

  dt <- fread(f)[score > 0] # turn the files into a data frame

# make sure the file is not empty

  if (nrow(dt) == 0) {
    message("Skipping ", fname, ": file is empty")
    next
  }

# ensure that the files being read have the seed_kmer column

  if (!"seed_kmer" %in% colnames(dt)) {
    message("Skipping ", fname, ": no seed_kmer column found")
    next
  }

# change all of the sequences into uppercase and grab the length of those sequences

  dna_vec <- toupper(as.character(dt$seed_kmer))
  lens    <- nchar(dna_vec)

# the calid sequences are greater than 1 in length and also not NaN
  valid_idx <- which(!is.na(dna_vec) & lens > 0)

  if (length(valid_idx) == 0) {
    message("Skipping ", fname, ": no valid sequences found")
    next
  }

# keep only the valid DNA sequences and their lengths
  dna_vec <- dna_vec[valid_idx]
  lens    <- lens[valid_idx]

# create a temporary data frame and calculate the instances of each nucleotide divided by the length of that sequence
  temp_dt <- data.table(
    A = str_count(dna_vec, "A") / lens,
    C = str_count(dna_vec, "C") / lens,
    G = str_count(dna_vec, "G") / lens,
    T = str_count(dna_vec, "T") / lens,
    sample = s_name,
    k_mer  = k_val
  )

# create a list of tables that will later go through an rbind
  all_data_list[[length(all_data_list) + 1]] <- temp_dt

# clean up the memory
  rm(dt, dna_vec, lens, temp_dt)
  gc()
}

# ensure that the table indeed has data and nothing went wrong in garbage collection
if (length(all_data_list) == 0) {
  stop("No valid data found after processing. Nothing to save or plot.")
}

# Rbind the tables to get master data frame
master_dt <- rbindlist(all_data_list, use.names = TRUE, fill = TRUE)

# progress statements
message("Final master_dt dimensions: ", nrow(master_dt), " rows x ", ncol(master_dt), " columns")
print(head(master_dt))
print(colnames(master_dt))

# Save
master_file <- file.path(output_folder, "kmer_GC_content_master_file.csv")
fwrite(master_dt, master_file)

message("Master table saved to: ", master_file)
message("Master table exists: ", file.exists(master_file))

# reorganize the data frame by grouping data together by ID and kmer length
plot_dt <- melt(
  master_dt,
  id.vars = c("sample", "k_mer"),
  measure.vars = c("A", "C", "G", "T"),
  variable.name = "Nucleotide",
  value.name = "Frequency"
)

message("Plot table dimensions: ", nrow(plot_dt), " rows x ", ncol(plot_dt), " columns")
print(head(plot_dt))

# make the box plot
p1 <- ggplot(plot_dt, aes(x = Nucleotide, y = Frequency, fill = Nucleotide)) +
  geom_boxplot(outlier.size = 0.5, color = "black", alpha = 0.8) +
  facet_grid(sample ~ k_mer) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_rect(fill = "gray90", color = "black"),
    panel.spacing = grid::unit(0.5, "lines")
  ) +
  labs(
    title = "Contigs Nucleotide Frequency Distribution",
    x = "Nucleotide",
    y = "Proportion (0-1)"
  )

# save
plot_file <- file.path(output_folder, "kmer_GC_content_representation.png")

png(filename = plot_file, width = 14, height = 10, units = "in", res = 300)
print(p1)
dev.off()

message("Plot saved to: ", plot_file)
message("Plot exists: ", file.exists(plot_file))

message("Success! Master CSV, PNG plot, are in: ", output_folder)

list.files(output_folder)