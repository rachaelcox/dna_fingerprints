#!/usr/bin/env Rscript

# 1. Load Libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(scales) # Required for labels = scales::percent

# 2. Ensure Output Directory Exists
# This is usually why ggsave fails on a server
if (!dir.exists("figures/images")) {
  dir.create("figures/images", recursive = TRUE)
}

# 3. Load and Clean
if (!file.exists("data/kmer_space/kmer_unique_report.csv")) {
  stop("Input file not found! Check your working directory.")
}

df_raw <- read.csv("data/kmer_space/kmer_unique_report.csv")

df_clean <- df_raw %>%
  mutate(
    Prop_Num = as.numeric(sub("%", "", Prop)) / 100,
    Cond = factor(Cond, levels = c("bg", "0w", "1w", "3w"))
  )

# 4. Create subsets
df_oligos <- df_clean %>% 
  filter(str_detect(Substrate, regex("oligo", ignore_case = TRUE)))

df_others <- df_clean %>% 
  filter(!str_detect(Substrate, regex("oligo", ignore_case = TRUE)))

# Shared Palette (Named to prevent 'No shared levels' warning)
pal <- c(
  "bg" = "grey30",
  "0w" = "#ffc1da",
  "1w" = "#e75480",
  "3w" = "#c90076",
  "#7a1fa2", "#1f8f6b", "#6ea8fe", "#4c72b0", "#353540","#DDA0DD"
)

# 5. Helper function
make_plot <- function(data_subset, plot_title) {
  if (nrow(data_subset) == 0) return(NULL) # Safety check for empty subsets
  
  plot_data <- data_subset %>%
    group_by(Substrate, Cond, K) %>%
    summarize(
      Mean_Prop = mean(Prop_Num, na.rm = TRUE), 
      SD_Prop = sd(Prop_Num, na.rm = TRUE), 
      .groups = 'drop'
    )
  
  ggplot(plot_data, aes(x = factor(K), y = Mean_Prop, fill = Cond)) +
    geom_col(position = position_dodge(0.9)) + 
    geom_errorbar(
      aes(ymin = Mean_Prop - SD_Prop, ymax = Mean_Prop + SD_Prop),
      position = position_dodge(0.9), width = 0.25
    ) +
    facet_wrap(~Substrate) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    scale_fill_manual(values = pal) +
    labs(title = plot_title, x = "K-mer Length (k)", 
         y = "Mean Proportion of Theoretical Space", fill = "Condition")
}

# 6. Generate and Save
# --- Plot All ---
p_all <- make_plot(df_clean, "Average K-mer Space Saturation (All)")
if(!is.null(p_all)) ggsave("figures/images/all_kmer_representation.png", p_all, width = 12, height = 8)

# --- Plot Oligos ---
p_oligos <- make_plot(df_oligos, "Average K-mer Space Saturation (Oligos)")
if(!is.null(p_oligos)) ggsave("figures/images/oligos_kmer_representation.png", p_oligos, width = 10, height = 6)

# --- Plot Others ---
p_others <- make_plot(df_others, "Average K-mer Space Saturation (Substrates)")
if(!is.null(p_others)) ggsave("figures/images/substrates_kmer_representation.png", p_others, width = 10, height = 6)

print("Script complete. Plots saved in figures/images/")