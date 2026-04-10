library(readr)
library(tidyverse)
pdf(NULL)

#   Figure 4 : Oligos Upset Plot                  
# ____________________________________________

library(ComplexUpset)

# read in the data after file combination
oligo1 <- read_csv("figures/combined/oligo1_data.csv")
oligo2 <- read_csv("figures/combined/oligo2_data.csv")
oligo3 <- read_csv("figures/combined/oligo3_data.csv")
oligo4 <- read_csv("figures/combined/oligo4_data.csv")

# mutate edges to substrates
fmt_data <- function(df, substrate){
  df <- df %>%
    select(edge) %>%
    mutate(sample = substrate) %>%
    mutate(k = nchar(edge)) %>%
    filter(k == 10)
  return(df)
}

# combine into one large df
all_edges <- bind_rows(
  oligo1 %>% fmt_data(substrate="Oligo 1"),
  oligo2 %>% fmt_data(substrate="Oligo 2"),
  oligo3 %>% fmt_data(substrate="Oligo 3"),
  oligo4 %>% fmt_data(substrate="Oligo 4")
)

# deduplicate edges and make df wider
edge_matrix <- all_edges %>%
  distinct(edge, sample) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

# install.packages("UpSetR")
library(UpSetR)

# define the set columns for upset plot
set_cols <- c("Oligo 1", "Oligo 2", "Oligo 3", "Oligo 4")

m <- as.data.frame(edge_matrix[, set_cols])

op <- par(col = "black")  # force text back to black

# Upset plot creation
oligos_upset <- UpSetR::upset(
  m,
  sets = set_cols,                
  keep.order = TRUE,              
  nsets = length(set_cols),
  nintersects = 20,
  order.by = "freq",
  
  main.bar.color = "#c90076",
  sets.bar.color = "grey80",
  matrix.colo =  "#c90076",
  point.size = 2.7,
  line.size = 0.7
)

# save figure
png(filename = "figures/images/oligos_upset.png", 
    width = 10, 
    height = 7, 
    units = "in", 
    res = 300)

print(oligos_upset)
dev.off()

# Figure 4 : Log2fc of kmer enrichment
#___________________________________________

# this is hardcoded in an would need to be changed when using something other than these data

oligo1 <- read_csv("figures/combined/oligo1_data.csv")
oligo2 <- read_csv("figures/combined/oligo2_data.csv")
oligo3 <- read_csv("figures/combined/oligo3_data.csv")
oligo4 <- read_csv("figures/combined/oligo4_data.csv")

# mutate labels to match the ouput of the label_kmers script
oligo_1_data_clean <- oligo1 %>%
  mutate(label = ifelse(label == "Perfect complement", "complement oligo match", "other"))

oligo_2_data_clean <- oligo2 %>%
  mutate(label = ifelse(label == "Perfect complement", "complement oligo match", "other"))

oligo_3_data_clean <- oligo3 %>%
  mutate(label = ifelse(label == "Perfect complement", "complement oligo match", "other"))

oligo_4_data_clean <- oligo4 %>%
  mutate(label = ifelse(label == "Perfect complement", "complement oligo match", "other"))

# hardcoding in the titles
oligo_1_data_clean$title <- "Oligo 1"
oligo_2_data_clean$title <- "Oligo 2"
oligo_3_data_clean$title <- "Oligo 3"
oligo_4_data_clean$title <- "Oligo 4"

# combining all of the data into one large df
combined_data <- rbind(
  oligo_1_data_clean,
  oligo_2_data_clean,
  oligo_3_data_clean,
  oligo_4_data_clean
)

# creation of log2fc plot, calling substrates the Timepoint
# Degust calculates the log2fc, this is then averaged across washes and substrates here
line_df <- combined_data %>%
  select(exp, kmer_size, label, `0w`, `1w`, `3w`) %>% 
  pivot_longer(
    cols = c(`0w`, `1w`, `3w`),
    names_to = "Timepoint",
    values_to = "log2FC"
  ) %>%
  group_by(exp, kmer_size, label, Timepoint) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE), # summary statistics here
    sd_log2FC = sd(log2FC, na.rm = TRUE),
    n = n(),
    se_log2FC = sd_log2FC / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(label_time = paste(label, Timepoint, sep = "_")) %>%
  filter(kmer_size <= 12)

lo2fc_figure <- ggplot(line_df, aes(x = as.factor(kmer_size), y = mean_log2FC,
                                    color = label_time, group = label_time)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_line(size = 1.5) + 
  geom_point(size = 2) +  # keep points
  geom_errorbar(
    aes(
      ymin = mean_log2FC - se_log2FC, # error bar creation here
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
      "other_0w" = "grey90",
      "other_1w" = "grey50",
      "other_3w" = "grey30"
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
    scales = "fixed"
  )

w = 10; h = 5
lo2fc_figure %>% ggsave("figures/images/lo2fc_figure.png", ., device = "png", 
                        width = w, height = h, units = "in")




# Figure 7 : Combined Substrate Upset Plot        
#______________________________________________

library(ComplexUpset)

oligo1 <- read_csv("figures/combined/oligo1_data.csv")
oligo2 <- read_csv("figures/combined/oligo2_data.csv")
oligo3 <- read_csv("figures/combined/oligo3_data.csv")
oligo4 <- read_csv("figures/combined/oligo4_data.csv")
streptavidin <- read_csv("figures/combined/streptavidin_data.csv")
glass <- read_csv("figures/combined/glass161_data.csv")
chitin <- read_csv("figures/combined/Chitin_Only_data.csv")


# function for modifying input data
fmt_data <- function(df, substrate){
  df <- df %>%
    select(edge) %>%
    mutate(sample = substrate) %>%
    mutate(k = nchar(edge)) %>%
    filter(k == 10)
  return(df)
}

all_edges <- bind_rows(
  oligo1 %>% fmt_data(substrate="Oligo 1"),
  oligo2 %>% fmt_data(substrate="Oligo 2"),
  oligo3 %>% fmt_data(substrate="Oligo 3"),
  oligo4 %>% fmt_data(substrate="Oligo 4"),
  streptavidin %>% fmt_data(substrate="Streptavidin"),
  glass %>% fmt_data(substrate="Glass"),
  chitin %>% fmt_data(substrate="Chitin")
)


edge_matrix <- all_edges %>%
  distinct(edge, sample) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)

# install.packages("UpSetR")
library(UpSetR)

# edge_matrix example columns: edge, Oligo 1, Oligo 2, Oligo 3, Oligo 4
set_cols <- c("Oligo 1", "Oligo 2", "Oligo 3", "Oligo 4", "Streptavidin", "Glass", "Chitin")

m <- as.data.frame(edge_matrix[, set_cols])

op <- par(col = "black")  # force text back to black

substrates_upset <- UpSetR::upset(
  m,
  sets = set_cols,                
  keep.order = TRUE,              
  nsets = length(set_cols),
  nintersects = 20,
  order.by = "freq",
  
  main.bar.color = "plum3",
  sets.bar.color = "grey80",
  matrix.color = "plum3",
  point.size = 2.7,
  line.size = 0.7
)

# save figure
png(filename = "figures/images/substrates_upset.png", 
    width = 10, 
    height = 7, 
    units = "in", 
    res = 300)

print(substrates_upset)
dev.off()


# Figure 8: NMDS Plot            
# _______________________________________________


library(vegan)
library(ggplot2)
library(patchwork)
library(tidyverse)

theme_set(theme_bw(base_size = 12))
options(scipen=100000)
pal <- c(
  "pink3",  # bubblegum pink
  "#c90076",  # raspberry pink
  "#b266ff",  # lavender purple
  "#7a1fa2",  # royal purple
  "#1f8f6b",  # deep jade green
  "#6ea8fe",  # periwinkle blue
  "#4c72b0",  # denim blue
  "#353540"   # charcoal
)

# reading in the 10mers cleaned file for plotting the NMDS
data_dir = "results"
exps <- c("oligo1", "oligo2", "oligo3", "oligo4",
          "chitin", "streptavidin", "glass")


df <- data.table::fread(file = "figures/clean_tables/all_substrates_enriched_10mer.csv")


df %>%
  mutate(na_count = rowSums(is.na(.))) %>%
  select(experiment, na_count) %>%
  filter(na_count > 0) 



set.seed(13)
df.dist <- vegdist(as.matrix(df %>% select(-experiment)), method = "bray") # creating the distance matrix using bray curtis
nmds <- metaMDS(df.dist, distance = "bray", k = 2, trymax = 100) # calculating NMDS ordinations for NMDS1 and 2, tries 100 random starting positions and chooses lowest stress

data.scores = as.data.frame(scores(nmds))
data.scores$experiment = df$experiment

# Calculating where the clusters are placed on the nmds plot, separated by substrate and wash

hull <- data.scores %>%
  extract(experiment, into = c("substrate", "wash"), regex = "(.*)_(.*)", remove = FALSE) %>% 
  group_by(substrate) %>%
  slice(chull(NMDS1, NMDS2))

pt_size = 3
p1 <- data.scores %>%
  extract(experiment, into = c("substrate", "wash"), regex = "(.*)_(.*)", remove = FALSE) %>%  # split on the underscore
  ggplot(aes(x = NMDS1, y = NMDS2, color = substrate)) +
  #geom_polygon(data = hull, alpha = 0.25, size = 1) +
  geom_point(size = pt_size+1, color = "black", aes(shape = wash)) +
  geom_point(size = pt_size, aes(shape = wash)) +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  annotate("text", 
           label = str_interp("stress = ${round(nmds$stress, 3)}"),
           x = -Inf, y = Inf,   
           hjust = -0.1,        
           vjust = 1.5)         


# Figure 9: Stress plot          
# ______________________________

s_result <- stressplot(nmds, plot = FALSE)

stress_df <- tibble(x = s_result$x, # true distances
                    y = s_result$y,
                    yf = s_result$yf) %>%
  pivot_longer(cols = c(y, yf),
               names_to = "var")


ord_dist <- stress_df %>% filter(var == "y") %>% pull(value) # pull out the rank ordered distances
fit_dist <- stress_df %>% filter(var == "yf") %>% pull(value) # pull out the distances predicted by the NMDS 2D representation

r2_linear <- cor(ord_dist, fit_dist)^2 # calculate linear R squared values
r2_nonmetric <- 1 - (sum((ord_dist - fit_dist)^2) / sum(ord_dist^2)) # calculate nonmetric fit R squared values


p2 <- stress_df %>%
  ggplot(aes(x = x, y = value)) +
  geom_point(data = stress_df %>% filter(var == "y"), alpha=0.5, size=2.5) +
  geom_step(data = stress_df %>% filter(var == "yf"),
            col = "#c90076", direction = "vh", size = 1) +
  labs(x = "Bray-Curtis dissimilarity", y = "Ordination distance") +
  
  annotate("text", x = min(stress_df$x)*0.9, y = max(stress_df$value), 
           label = paste0("Non-metric fit, R^2 = ", round(r2_nonmetric, 3)), 
           hjust = 0) +
  annotate("text", x = min(stress_df$x)*0.9, y = max(stress_df$value)*0.9, 
           label = paste0("Linear fit, R^2 = ", round(r2_linear, 3)), 
           hjust = 0)


NMDS_stressplot <- p1 + p2

w = 10; h = 5
NMDS_stressplot %>% ggsave("figures/images/NMDS_stressplot.png", ., device = "png", 
                           width = w, height = h, units = "in")

# Figure 8b PERMANOVA
#______________________

# Prepare the meta data by splittign on substrate and wash
df_meta <- df %>%
  select(experiment) %>%
  separate(experiment, into = c("substrate", "wash"), sep = "_", remove = FALSE)

# run permanova by substrate and shuffle data for Null hypothesis 999 times
permanova_substrate <- adonis2(df.dist ~ substrate, data = df_meta, permutations = 999)

# creates the dispersion objects to determine within group dispersion and calculate the ANOVA of the group dispersion
disp_substrate <- betadisper(df.dist, df_meta$substrate)
disp_anova <- anova(disp_substrate)

# calculate ANOSIM in the same way - to tell use within group 
anosim_substrate <- anosim(df.dist, grouping = df_meta$substrate, permutations = 999)

print("PERMANOVA by Substrate:")
print(permanova_substrate)
print("--- BETA DISPERSION RESULTS ---")
print(disp_anova)
print("ANOSIM by Substrate:")
summary(anosim_substrate)


capture.output({

  cat("===== PERMANOVA by Substrate =====\n")
  print(permanova_substrate)

  cat("\n===== BETA DISPERSION RESULTS =====\n")
  print(disp_anova)

  cat("\n===== ANOSIM by Substrate =====\n")
  print(summary(anosim_substrate))

}, file = "figures/combined/NMDS_statistics.txt")