# SPACE APTAMERS FIGURE GENERATION

# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ 
# ♡              libraries                ♡ 
# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ 

library(tidyverse)
library(patchwork)
library(fs)
library(ggVennDiagram)
library(ggConvexHull)
library(ComplexUpset)
library(dplyr)
library(tidyr)
library(ggplot2)

# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ 
# ♡               data files              ♡ 
# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ 

strep_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/strep.csv")
oligo_1_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/oligo_1.csv")
oligo_2_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/oligo_2.csv")
oligo_3_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/oligo_3.csv")
oligo_4_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/oligo_4.csv")
Chitin_Only_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/chitin.csv")
glass161_data <- read_csv("/stor/work/Marcotte/project/zoya/oligos/figure_generation/glass.csv")
k_space_file = '/stor/work/Marcotte/project/rmcox/oligos/results/oligo_kmer_space_sampled.csv'
oligo1 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_1_data.csv")
oligo2 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_2_data.csv")
oligo3 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_3_data.csv")
oligo4 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_4_data.csv")


pinkpal = c("#353540", "mistyrose", "lightpink", "#C90076")


df <- read_csv(k_space_file) %>%
  mutate(total_possible_kmers = 4^kmer_size,
         #k_label = glue::glue("{kmer_size}\n({total_possible_kmers})")
         k_label = glue::glue("{kmer_size}\n({formatC(total_possible_kmers, format='e', digits=1)})")
  )

df$k_label <- fct_reorder(df$k_label, df$kmer_size)
p1 <- df %>%
  ggplot(aes(x=k_label, y=proportion_sampled, fill=as.factor(n_washes))) +
  geom_col(position = position_dodge2(preserve='single', width=0.21)) +
  facet_wrap(~exp, ncol=1) +
  scale_fill_manual(values=pinkpal) +
  labs(x = "Size of k (total possible unique sequences of k)",
       y = "Proportion of sequence space covered",
       fill = "# of washes") +
  theme(legend.position = "top")

p1

# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡  
# ♡  Figure 5 : Upset Plot                  ♡ 
# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡  

# function for modifying input data
fmt_data <- function(df, substrate){
  df <- df %>%
    select(edge) %>%
    mutate(sample = substrate) %>%
    mutate(k = nchar(edge)) #%>%
  # filter(k == 10)
  return(df)
}

all_edges <- bind_rows(
  oligo1 %>% fmt_data(substrate="Oligo 1"),
  oligo2 %>% fmt_data(substrate="Oligo 2"),
  oligo3 %>% fmt_data(substrate="Oligo 3"),
  oligo4 %>% fmt_data(substrate="Oligo 4")
)


edge_matrix <- all_edges %>%
  distinct(edge, sample) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

p3 <- ComplexUpset::upset(
  edge_matrix,
  intersect = c(
    "Oligo 1", "Oligo 2", "Oligo 3", "Oligo 4"
  ),
  base_annotations=list(
    'Intersection size'=intersection_size(
      fill="#c90076",
      color="#353540",
      text=list(
        size=3,
        #vjust=-.01,
        hjust=-0.0005,
        angle=45
      )
    )
  ),
  n_intersections=20,
  width_ratio=0.15
)
p3

# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡  
# ♡  Figure 6 : Log2fc of kmer enrichment   ♡ 
# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡

oligo_1_data_clean <- oligo_1_data %>%
  mutate(label = ifelse(label == "complement oligo match", "complement oligo match", "other"))

oligo_2_data_clean <- oligo_2_data %>%
  mutate(label = ifelse(label == "complement oligo match", "complement oligo match", "other"))

oligo_3_data_clean <- oligo_3_data %>%
  mutate(label = ifelse(label == "complement oligo match", "complement oligo match", "other"))

oligo_4_data_clean <- oligo_4_data %>%
  mutate(label = ifelse(label == "complement oligo match", "complement oligo match", "other"))

oligo_1_data_clean$title <- "Oligo 1"
oligo_2_data_clean$title <- "Oligo 2"
oligo_3_data_clean$title <- "Oligo 3"
oligo_4_data_clean$title <- "Oligo 4"


combined_data <- rbind(
  oligo_1_data_clean,
  oligo_2_data_clean,
  oligo_3_data_clean,
  oligo_4_data_clean
)

line_df <- combined_data %>%
  select(exp, kmer_size, label, `0w`, `1w`, `3w`) %>% 
  pivot_longer(
    cols = c(`0w`, `1w`, `3w`),
    names_to = "Timepoint",
    values_to = "log2FC"
  ) %>%
  group_by(exp, kmer_size, label, Timepoint) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    sd_log2FC = sd(log2FC, na.rm = TRUE),
    n = n(),
    se_log2FC = sd_log2FC / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(label_time = paste(label, Timepoint, sep = "_")) %>%
  filter(kmer_size <= 12)

p4 <- ggplot(line_df, aes(x = as.factor(kmer_size), y = mean_log2FC,
                          color = label_time, group = label_time)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_line(size = 1.5) + 
  geom_point(size = 2) +  # keep points
  geom_errorbar(
    aes(
      ymin = mean_log2FC - se_log2FC,
      ymax = mean_log2FC + se_log2FC
    ),
    width = 0.3,
    color = "black",
    linetype = "solid"
  ) +
  scale_color_manual(
    values = c(
      "complement oligo match_0w" = "#ffc1da",
      "complement oligo match_1w" = "#e75480",
      "complement oligo match_3w" = "#c90076",
      "other_0w" = "grey80",
      "other_1w" = "grey40",
      "other_3w" = "black"
    ),
    labels = c(
      "complement oligo match_0w" = "Complementary to target oligo (0w)",
      "complement oligo match_1w" = "Complementary to target oligo (1w)",
      "complement oligo match_3w" = "Complementary to target oligo (3w)",
      "other_0w" = "All other K-mer sequences (0w)",
      "other_1w" = "All other K-mer sequences (1w)",
      "other_3w" = "All other K-mer sequences (3w)"
    )
  ) +
  labs(
    #title = "Mean K-mer log2 Fold Change by Wash Condition",
    x = "K-mer size",
    y = "Log2fc (mean)",
    color = "Sequence category"
  ) +
  #theme_minimal(base_size = 16) +
  theme_bw() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold")
  ) +
  facet_wrap(
    ~ exp,
    nrow = 1,
    scales = "fixed",
    labeller = as_labeller(c(
      "oligo 1" = "Oligo 1",
      "oligo 2" = "Oligo 2",
      "oligo 3" = "Oligo 3",
      "oligo b" = "Oligo 4"
    ))
  )
print(p4)



# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡  
# ♡  Final Figure: Combined panels          ♡ 
# ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡ ♡
theme_set(theme_bw(base_size = 12))
options(scipen=100000)

finale <- (p1 + (p3 / p4 + plot_layout(heights = c(1, 1)))) +
  plot_layout(widths = c(0.75, 1.25), guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(
    axis.title = element_text(size = 8),
    legend.position = "bottom"
  )


print(finale)

ggsave("spapt_fig4.pdf", plot = finale, width = 20, height = 20, units = "in")


