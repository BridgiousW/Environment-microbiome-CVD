#######
library(tidyverse)
library(igraph)
library(ggraph)

# Step 1: Create the dataset
data <- tribble(
  ~log2FoldChange, ~padj, ~Abundant_Group, ~taxa, ~Var2, ~value, ~absCorr, ~associated_CVD_risk_factor,
  2.1867, 2.58e-12, "Positive", "Altererythrobacter", "HMDB31050", -0.829, 0.829, "DIASTOLIC BP",
  2.1129, 7.28e-11, "Positive", "Arthrobacter", "HMDB31050", -0.847, 0.847, "DIASTOLIC BP",
  1.3364, 9.50e-06, "Positive", "Devosia", "HMDB31050", -0.819, 0.819, "DIASTOLIC BP",
  1.4811, 6.24e-07, "Positive", "Domibacillus", "HMDB32627", -0.800, 0.800, "DIASTOLIC BP",
  2.1777, 9.73e-11, "Positive", "Ellin6055", "HMDB31050", -0.821, 0.821, "DIASTOLIC BP",
  1.0694, 0.00062, "Positive", "Geodermatophilus", "HMDB31050", -0.813, 0.813, "DIASTOLIC BP",
  1.0117, 0.00011, "Positive", "Kapabacteriales", "HMDB32627", -0.767, 0.767, "DIASTOLIC BP",
  1.5948, 1.83e-07, "Positive", "Kribbella", "HMDB31050", -0.797, 0.797, "DIASTOLIC BP",
  1.0321, 0.00061, "Positive", "Longimicrobiaceae", "HMDB31050", -0.798, 0.798, "DIASTOLIC BP",
  2.6592, 2.48e-16, "Positive", "Lysobacter", "HMDB31050", -0.853, 0.853, "DIASTOLIC BP",
  1.0511, 0.00337, "Positive", "Nitrospira", "HMDB31050", -0.810, 0.810, "DIASTOLIC BP",
  1.8145, 4.45e-07, "Positive", "Pseudarthrobacter", "HMDB31050", -0.825, 0.825, "DIASTOLIC BP",
  3.8699, 8.84e-22, "Positive", "Vicinamibacteraceae", "HMDB32627", -0.777, 0.777, "DIASTOLIC BP",
  2.2985, 2.80e-11, "Positive", "Gaiella", "0.95_764.7010n", -0.671, 0.671, "LDL CHOLESTEROL",
  2.2985, 2.80e-11, "Positive", "Gaiella", "HMDB56087", -0.655, 0.655, "LDL CHOLESTEROL",
  1.6086, 3.99e-08, "Positive", "Arenimonas", "0.95_764.7010n", -0.654, 0.654, "LDL CHOLESTEROL"
)

# Step 2: Create edges
edges_microbe_metab <- dplyr::select(data, from = taxa, to = Var2)
edges_metab_cvd <- dplyr::select(data, from = Var2, to = associated_CVD_risk_factor)
edges <- bind_rows(edges_microbe_metab, edges_metab_cvd)

# Step 3: Create node types
nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(type = case_when(
    name %in% data$taxa ~ "Microbe",
    name %in% data$Var2 ~ "Metabolite",
    name == "LDL CHOLESTEROL" ~ "LDL",
    name == "DIASTOLIC BP" ~ "DIAS"
  ))

# Step 4: Set colors
node_colors <- c(
  "Microbe" = "#08306B",        
  "Metabolite" = "#FFD700",     
  "LDL" = "#228B22",            
  "DIAS" = "#B22222"            
)

# Step 5: Create and plot graph
g <- graph_from_data_frame(edges, vertices = nodes)

ggraph(g, layout = "fr") +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm'),
                 edge_colour = "gray60") +
  geom_node_point(aes(color = type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
  scale_color_manual(values = node_colors, name = "Node Type",
                     labels = c("Microbe", "Metabolite", "LDL CHOLESTEROL", "DIASTOLIC BP")) +
  theme_void() +
  labs(title = "Microbe–Metabolite–CVD Risk Factor Network") +
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold"))


