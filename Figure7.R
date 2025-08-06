#####
library(tidyverse)
library(igraph)
library(ggraph)
setwd("/Users/mbmhvbw2/metabolomics")
merged_all<-read.csv("merged_all.csv")
# Start with your dataset
data <- merged_all  # Replace with your actual dataframe name if different

# Add correlation direction based on the `value` column (assuming that's the signed correlation)
edges_microbe_metab <- data %>%
  dplyr::select(from = taxa, to = Var2, corr = value) %>%  # 'value' is signed correlation
  mutate(corr_dir = ifelse(corr > 0, "Positive", "Negative"))

# Create metabolite → CVD marker edges (no correlation direction here)
edges_metab_cvd <- data %>%
  dplyr::select(from = Var2, to = marker) %>%
  mutate(corr_dir = NA)
##
###
# Load required libraries
library(tidyverse)
library(igraph)
library(ggraph)
library(scales)

# Load your dataset (assuming it's named `data`)
# If not already loaded:
# data <- read.csv("your_file.csv")  # replace with actual file if needed

# Step 1: Prepare edge lists
edges_microbe_metab <- data %>%
  select(from = taxa, to = Var2)

edges_metab_cvd <- data %>%
  select(from = Var2, to = marker)

# Combine edges
edges <- bind_rows(edges_microbe_metab, edges_metab_cvd)

# Step 2: Determine correlation direction (positive/negative) for microbe-metabolite edges
edges <- edges %>%
  left_join(data %>% dplyr::select(taxa, Var2, value), by = c("from" = "taxa", "to" = "Var2")) %>%
  mutate(corr_dir = case_when(
    !is.na(value) & value > 0 ~ "Positive",
    !is.na(value) & value < 0 ~ "Negative",
    TRUE ~ "None"
  ))

# Step 3: Create node metadata
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(type = case_when(
    name %in% data$taxa ~ "Microbe",
    name %in% data$Var2 ~ "Metabolite",
    name %in% data$marker ~ "CVD Risk"
  ))

# Step 4: Create igraph object
g <- graph_from_data_frame(edges, vertices = nodes)

# Step 5: Plot network using ggraph
ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = corr_dir),
                 arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_width = 0.8,
                 show.legend = TRUE,
                 na.rm = TRUE) +
  geom_node_point(aes(color = type), size = 6, na.rm = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, na.rm = TRUE) +
  scale_color_manual(values = c("Microbe" = "skyblue",
                                "Metabolite" = "orange",
                                "CVD Risk" = "lightgreen")) +
  scale_edge_color_manual(values = c("Positive" = "blue",
                                     "Negative" = "red",
                                     "None" = "grey70")) +
  theme_void() +
  ggtitle("Microbe–Metabolite–CVD Risk Factor Network") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"))


###