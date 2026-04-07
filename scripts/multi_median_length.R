# setwd to oligos_redux (or the main data directory)

library(Biostrings)
library(ggplot2)
library(stringr)
library(dplyr)

# create color palette for 5-15 kmers
pal <- c(
    "pink3",# bubblegum pink
    "#C90076",# raspberry pink
    "#B266FF",# lavender purple  
    "#7A1FA2",# royal purple 
    "#05cd94", # light green
    "#1F8F6B",# deep jade green  
    "#04452fe5", # darker green
    "#12341e", # darkgreen
    "#6EA8FE",# periwinkle blue  
    "#4C72B0",# denim blue  
    "#444476"  # darker blue 
)

# Path to your processed data folder
parent_folder <- "results"

# Find all FASTA files recursively across all subfolders
file_list <- list.files(path = parent_folder, 
                        pattern = "\\.(fasta|fa)$", 
                        full.names = TRUE, 
                        recursive = TRUE)

# Process files and extract the "Sample Type" from two directory levels up named for the substrate
all_data <- do.call(rbind, lapply(file_list, function(f) {
  # Read sequences
  dna <- readDNAStringSet(f)
  sample_name <- basename(dirname(dirname(f)))
  
  # Extract KmerSize from filename before 'mers'
  kmer_size <- as.numeric(str_extract(basename(f), "\\d+(?=mers)"))
  

  # calculate the length of the dna
  data.frame(
    SampleType = sample_name,
    KmerSize = kmer_size,
    SeqLength = width(dna)
  )
}))

# order the x axis as ascending when creating the plot
all_data <- all_data %>%
  mutate(PlotLabel = factor(KmerSize, levels = sort(unique(KmerSize))))


# making the plots

# Get a list of every unique substrate
unique_samples <- unique(all_data$SampleType)

# Loop through each substrate to create separate boxplots and CSVs of median sequence length
for (s in unique_samples) {
  
  #  create a sample subset for each substrate, so this plots each substrate individually
  sample_subset <- all_data %>% filter(SampleType == s)
  
  #  Create the box plot
  p <- ggplot(sample_subset, aes(x = PlotLabel, y = SeqLength, fill = PlotLabel)) +
    geom_boxplot(color = "black") + 
    scale_fill_manual(values = pal) + 
    labs(
      title = paste("Contiguous Sequence Length Distribution:", s), # Correct Label
      x = "k-mer length",
      y = "Sequence Length (bp)"
    ) +
    theme_bw() + 
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  #  Save the Plot
  chart_filename <- paste0(s, "_fasta_chart_boxplots.png")
  # create.dir = TRUE ensures R makes the folder if it's missing
  ggsave(filename = chart_filename, plot = p, path = "figures/images", 
         create.dir = TRUE, width = 8, height = 5, dpi = 300)
  
  #  Calculate and Save the Median CSV for this specific sample
  sample_summary <- sample_subset %>%
    group_by(KmerSize) %>%
    summarize(MedianLength = median(SeqLength, na.rm = TRUE)) %>%
    arrange(desc(MedianLength))
  
  csv_filename <- paste0(s, "_kmer_median_lengths.csv")
  write.csv(sample_summary, file = file.path("consensus/median_lengths", csv_filename), row.names = FALSE)
  
  print(paste("Finished processing:", s))
}