library(vegan)
library(ggplot2)
library(patchwork)
library(ComplexUpset)

theme_set(theme_bw(base_size = 12))
options(scipen=100000)
pal <- c( "#8A2BE2", "#2986cc", "#00A201", 
          "#c90076", "#E67E22", "#EDAE49",
          "#353540", "#FB3640")


# upset plot -------------------------------------------------------------------

# read in data
chitin <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/Chitin_Only_data.csv")
glass <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/glass161_data.csv")
strep <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/strep_data.csv")
oligo1 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_1_data.csv")
oligo2 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_2_data.csv")
oligo3 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_3_data.csv")
oligo4 <- read_csv("/stor/work/Ellington/sanchita/oligos/results/pca/oligo_4_data.csv")

# function for modifying input data
fmt_data <- function(df, substrate){
  df <- df %>%
    select(edge) %>%
    mutate(sample = substrate) %>%
    mutate(k = nchar(edge)) %>%
    filter(k == 10)
  return(df)
}

# combine all 'edge' values with sample labels
all_edges <- bind_rows(
  chitin %>% fmt_data(substrate="Chitin"),
  glass %>% fmt_data(substrate="Glass"),
  strep %>% fmt_data(substrate="Streptavidin"),
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
    "Chitin", "Glass", "Streptavidin", 
    "Oligo 1", "Oligo 2", "Oligo 3", "Oligo 4"
  ),
  base_annotations=list(
    'Intersection size'=intersection_size(
      fill="#8A2BE2",
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

# nmds plot --------------------------------------------------------------------

data_dir = "oligos/results"
exps <- c("oligo_1", "oligo_2", "oligo_3", "oligo_b",
          "chitin", "streptavidin", "glass")
df <- data.table::fread(file = "oligos/results/clean_tables/all_substrates_10mers.csv")

df %>%
  mutate(na_count = rowSums(is.na(.))) %>%
  select(experiment, na_count)

set.seed(13)
df.dist <- vegdist(as.matrix(df %>% select(-experiment)), method = "bray")
nmds <- metaMDS(df.dist, distance = "bray", k = 2, trymax = 100)

data.scores = as.data.frame(scores(nmds))
data.scores$experiment = df$experiment

hull <- data.scores %>%
  separate(experiment, into = c("substrate", "wash"), sep = "_") %>% 
  group_by(substrate) %>%
  slice(chull(NMDS1, NMDS2))

pt_size = 3
p1 <- data.scores %>%
  separate(experiment, into = c("substrate", "wash"), sep = "_") %>% 
  mutate(substrate = case_when(
    substrate == "oligo1" ~ "Oligo 1",
    substrate == "oligo2" ~ "Oligo 2",
    substrate == "oligo3" ~ "Oligo 3",
    substrate == "oligo4" ~ "Oligo 4",
    substrate == "glass" ~ "Glass",
    substrate == "strep" ~ "Streptavidin",
    substrate == "chitin" ~ "Chitin"
  )) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = substrate)) +
  #geom_polygon(data = hull, alpha = 0.25, size = 1) +
  geom_point(size = pt_size+1, color = "black", aes(shape = wash)) +
  geom_point(size = pt_size, aes(shape = wash)) +
  scale_color_manual(values=pal) +
  scale_fill_manual(values=pal) +
  annotate("text", label = str_interp("stress = ${round(nmds$stress, 3)}"),
           x = max(data.scores$NMDS1)*0.5, y = max(data.scores$NMDS2)) +
  theme(legend.title = element_blank())
p1

# stress plot (how to make pretty?)
stressplot(nmds, l.col = "#c90076", p.col = "#353540")

# shepherd plot ----------------------------------------------------------------
# this looks promising: https://stackoverflow.com/questions/69216955/ggplot2-version-of-shepard-plot-i-e-veganstressplot
str(stressplot(nmds))
# create a tibble that contains the data from stressplot
stress_df <- tibble(x = stressplot(nmds)$x,
                    y = stressplot(nmds)$y,
                    yf = stressplot(nmds)$yf) %>%
  # change data to long format
  pivot_longer(cols = c(y, yf),
               names_to = "var")

p2 <- stress_df %>%
  ggplot(aes(x = x, y = value)) +
  geom_point(data = stress_df %>% filter(var == "y"), alpha=0.5, size=2.5) +
  geom_step(data = stress_df %>% filter(var == "yf"),
            col = "#c90076", direction = "vh", size = 1) +
  labs(x = "Bray-Curtis dissimilarity", y = "Ordination distance") +
  annotate("text", label = "Non-metric fit, R^2 = 0.992", hjust = 0,
           x = min(stress_df$x)*0.9, y = max(stress_df$value)) +
  annotate("text", label = "Linear fit, R^2 = 0.966", hjust = 0,
           x = min(stress_df$x)*0.9, y = max(stress_df$value)*0.9)
p2

# full panel -------------------------------------------------------------------

panel <- (p3 / (p1 + p2)) + plot_annotation(tag_levels = 'A')
panel

w = 12; h = 9
panel %>% ggsave("oligos/figures/panel_all_substrates_10mers.png", ., device = "png",
                 width = w, height = h, units = "in")
