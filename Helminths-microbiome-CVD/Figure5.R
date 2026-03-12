
library(ggplot2)
library(ggrepel)

require(mt)
require(ggplot2)
library(ropls)
library(pROC)

##### Panel A #########
VolcanoplotdataPos_NEG <- read.csv("~/metabolomics/VolcanoplotdataPos_NEG.csv")

VolcanoplotdataPos_NEG$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
VolcanoplotdataPos_NEG$diffexpressed[VolcanoplotdataPos_NEG$Log2FC > 1.0 & VolcanoplotdataPos_NEG$p.value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
VolcanoplotdataPos_NEG$diffexpressed[VolcanoplotdataPos_NEG$Log2FC < -1.0 & VolcanoplotdataPos_NEG$p.value < 0.05] <- "DOWN"
mycolors <- c("#D55E00", "#0072B2", "#009E73")
names(mycolors) <- c("DOWN", "UP", "NO")

merged_df_selected <- readRDS('~/metabolomics/merged_df_selected.rds')

VolcanoplotdataPos_NEG$HMBD_ID<-merged_df_selected$chemical_ID[match(VolcanoplotdataPos_NEG$metabolite, merged_df_selected$Compound)]
VolcanoplotdataPos_NEG$HMBD_ID[is.na(VolcanoplotdataPos_NEG$HMBD_ID)] <- VolcanoplotdataPos_NEG$metabolite[is.na(VolcanoplotdataPos_NEG$HMBD_ID)]
##
# Top 10 UP
top_up <- VolcanoplotdataPos_NEG[VolcanoplotdataPos_NEG$diffexpressed == "UP", ]
top_up <- top_up[order(top_up$p.value), ][1:min(10, nrow(top_up)), ]

# Top 10 DOWN
top_down <- VolcanoplotdataPos_NEG[VolcanoplotdataPos_NEG$diffexpressed == "DOWN", ]
top_down <- top_down[order(top_down$p.value), ][1:min(10, nrow(top_down)), ]

# Combine
top_labels <- rbind(top_up, top_down)


p <- ggplot(data = VolcanoplotdataPos_NEG, 
            aes(x = Log2FC, y = -log10(p.value), color = diffexpressed)) +
  geom_point() +
  geom_text_repel(data = top_labels,
                  aes(label = HMBD_ID),
                  size = 3, max.overlaps = 100) +
  geom_vline(xintercept = c(-1.0, 1.0), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  scale_color_manual(values = mycolors) +
  theme_minimal() +
  labs(color = "Differential Abundance")

p

VolcanoplotdataPos_NEG$HMBD_ID<-merged_df_selected$chemical_ID[match(VolcanoplotdataPos_NEG$metabolite, merged_df_selected$Compound)]
VolcanoplotdataPos_NEG$HMBD_ID[is.na(VolcanoplotdataPos_NEG$HMBD_ID)] <- VolcanoplotdataPos_NEG$metabolite[is.na(VolcanoplotdataPos_NEG$HMBD_ID)]
##
# Top 10 UP
top_up <- VolcanoplotdataPos_NEG[VolcanoplotdataPos_NEG$diffexpressed == "UP", ]
top_up <- top_up[order(top_up$p.value), ][1:min(10, nrow(top_up)), ]

# Top 10 DOWN
top_down <- VolcanoplotdataPos_NEG[VolcanoplotdataPos_NEG$diffexpressed == "DOWN", ]
top_down <- top_down[order(top_down$p.value), ][1:min(10, nrow(top_down)), ]

# Combine
top_labels <- rbind(top_up, top_down)

mycolors <- c("DOWN" = "#D55E00", "UP" = "#0072B2", "NO" = "#009E73")

p <- ggplot(data = VolcanoplotdataPos_NEG, 
            aes(x = Log2FC, y = -log10(p.value), color = diffexpressed)) +
  geom_point() +
  geom_text_repel(data = top_labels,
                  aes(label = HMBD_ID),
                  size = 3, max.overlaps = 100) +
  geom_vline(xintercept = c(-1.0, 1.0), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  scale_color_manual(values = mycolors) +
  theme_minimal() + annotate("text", x = -1.5, y = 4, label = "Uninfected", size = 5) +
  annotate("text", x = 1.5, y = 4, label = "Infected", size = 5) +
  labs(color = "Differential Abundance")
p

##### Panel B #########

DataforPLSDA<- read.csv("~/metabolomics/DataforPLSDA.csv",check.names=FALSE)

DataforPLSDA$Hel <- meta_data$Hel_inf_status[
  match(DataforPLSDA$NAME, meta_data$slabno)
]
table(DataforPLSDA$Hel)

DataforPLSDA$case2 <- factor(DataforPLSDA$case, levels = c(0,1), labels = c('Uinfected', 'Infected') )
table(DataforPLSDA$case2)  

plsdamodel<-plslda(case2~., data = DataforPLSDA, pls="oscorespls")
#summary(plsdamodel)
summary(plsdamodel$pls.out)

comp1 <- plsdamodel$pls.out$scores[,1]
comp2 <- plsdamodel$pls.out$scores[,2]

df2 <- data.frame(comp1, comp2, case=DataforPLSDA$case2)

p <- ggplot(data = df2, aes(x=comp1, y=comp2, fill = case, color=case))+
  geom_point()+theme(axis.title.x = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  scale_color_manual("case",values=rev(c("#DAA520","#4682B4"))) +
  theme(axis.text.x = element_text(angle=0)
        ,panel.background = element_blank(),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size=15),#
        panel.border = element_rect(fill=NA),
        title = element_text(size = 12, face = "bold")
  )+stat_ellipse() +
  xlab("Component 1 (8.74%)")+ylab("Component 2 (7.76%)")
p


##### Panel C #########
# Load libraries
library(ggplot2)
library(dplyr)
library(viridis) # for the beautiful viridis color scale

# Create the dataframe
df <- data.frame(
  Pathway = c(
    "NR1H2 & NR1H3 regulate cholesterol uptake",
    "NR1H2 & NR1H3 regulate bile acid homeostasis",
    "NR1H2 & NR1H3 regulate gluconeogenesis",
    "NR1H2 & NR1H3 regulate lipogenesis",
    "NR1H2 & NR1H3 regulate triglyceride lipolysis",
    "NR1H2 & NR1H3 regulate cholesterol transport",
    "Cholesterol biosynthesis (skeletal dysplasias)"
  ),
  P_value = c(0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0003, 0.0011)
)

# Prepare data
df <- df %>%
  mutate(logP = -log10(P_value)) %>%
  arrange(P_value)

# Make the bar plot with color mapped to p-value
ggplot(df, aes(x = logP, y = reorder(Pathway, logP), fill = logP)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "rocket", direction = -1) +
  labs(
    x = expression(-log[10](P-value)),
    y = "Biological Pathway",
    title = "Biological pathways enriched by metabolites",
    fill = expression(-log[10](P-value))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


