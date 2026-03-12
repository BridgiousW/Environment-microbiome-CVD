require(phyloseq)
require(splitstackshape)
require(vegan)
require(ggplot2)
require(pheatmap)
require(DESeq2)
require(dunn.test)
require(cluster)
require(factoextra)
require(dplyr)
require(vegan)
require(gplots)
require(tidyr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(dplyr)

physeq<-readRDS('~/Novogene/codeNaturecomms/PhyloseqObjectHelminths.rds')
physeq_genus<-taxa_level(physeq, which_level='Genus') # see helper functions below for taxa_level function
otu_table(physeq_genus)<-otu_table(physeq_genus)+1
dds<-phyloseq_to_deseq2(physeq_genus, ~ Hel_inf_status)
dds<-DESeq(dds, test="Wald", fitType="local")
res1<-results(dds, cooksCutoff = FALSE)
alpha<-0.05
sigtab<-res1[which(res1$padj < alpha), ]
sigtab<-sigtab[which(abs(sigtab$log2FoldChange)>=1.0), ]
sigtab<-as.data.frame(sigtab)
sigtab$taxa <- rownames(sigtab)

## INSULIN
cp_cp_insulin_longl_long_microbes_filtered<-read.csv('~/Novogene/codeNaturecomms/cp_cp_insulin_longl_long_microbes_filtered.csv')
merged_data_insulin <- sigtab %>%
  inner_join(cp_cp_insulin_longl_long_microbes_filtered, by = c("taxa" = "Var1"))

# LDL
cp_LDLchol_long_microbes_filtered<-read.csv('~/Novogene/codeNaturecomms/cp_LDL_long_microbes_filtered.csv')
merged_data_LDL <- sigtab %>%
  inner_join(cp_LDLchol_long_microbes_filtered, by = c("taxa" = "Var1"))

# T_COL
cp_totalchol_long_microbes_filtered <- read.csv('~/Novogene/codeNaturecomms/cp_totalchol_long_microbes_filtered.csv')
merged_data_Totalchol <- sigtab %>%
  inner_join(cp_totalchol_long_microbes_filtered, by = c("taxa" = "Var1"))

# Sys
cp_cp_sys_long_microbes_filtered <- read.csv('~/Novogene/codeNaturecomms/cp_cp_sys_long_microbes_filtered.csv')
merged_data_systolic <- sigtab %>%
  inner_join(cp_cp_sys_long_microbes_filtered, by = c("taxa" = "Var1"))

cp_dias_long_microbes_filtered <- read.csv('~/Novogene/codeNaturecomms/cp_dias_long_microbes_filtered.csv')
merged_data_dias <- sigtab %>%
  inner_join(cp_dias_long_microbes_filtered, by = c("taxa" = "Var1"))


####
merged_all <- bind_rows(
  merged_data_LDL %>% mutate(marker = "LDL"),
  merged_data_Totalchol %>% mutate(marker = "Total Cholesterol"),
  merged_data_dias %>% mutate(marker = "Diastolic BP"),
  merged_data_systolic %>% mutate(marker = "Systolic BP"),
  merged_data_insulin %>% mutate(marker = "Insulin")
)

####
# and has columns: Microbe, Correlation (can be negative)
# For each Microbe, get the row with the highest absolute correlation
top_per_microbe <- merged_data_dias %>%
  group_by(taxa) %>%
  slice_max(order_by = absCorr , n = 1, with_ties = FALSE) %>%
  ungroup()
###
top_per_microbe$associated_CVD_risk_factor <- "DIASTOLIC BP"
merged_data_LDL$associated_CVD_risk_factor <- "LDL CHOLESTEROL"
finaldataSchistomicrobemetabolome <- rbind(top_per_microbe,merged_data_LDL )

# Add a binary indicator column for CVD-associated taxa
cvd_taxa <- c("Altererythrobacter", "Arthrobacter","Devosia","Domibacillus","Ellin6055","Geodermatophilus","Kapabacteriales","Kribbella","Longimicrobiaceae","Lysobacter","Nitrospira","Pseudarthrobacter","Vicinamibacteraceae","Gaiella","Arenimonas")  # <-- example list

df<-data.frame(sample_data(physeqcasegenus))["Hel_inf_status"]
sigtab$Abundant_Group <- levels(as.factor(df$Hel_inf_status))[as.numeric(sigtab$log2FoldChange>0) + 1]
df_wf<-sigtab[c("log2FoldChange","padj", "Abundant_Group")]
df_wf$taxa<-rownames(df_wf)
df_wf$taxa<-gsub("_"," ", df_wf$taxa)
df_wf<-df_wf[base::order(-df_wf$log2FoldChange),]

data <- df_wf %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = ifelse(padj < 0.05, "Significant", "Not Significant")
  )


data_clean <- data %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = ifelse(padj < 0.05, "Significant", "Not Significant"),
    cvd_associated = ifelse(taxa %in% cvd_taxa, "Yes", "No")
  )

# Volcano plot with shape mapped to CVD association
ggplot(data_clean, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = log2FoldChange, shape = cvd_associated), size = 3) +
  geom_text_repel(
    data = filter(data_clean, significance == "Significant"),
    aes(label = taxa),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
  scale_shape_manual(values = c("Yes" = 17, "No" = 16)) +  # triangle for CVD, circle for others
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Differential Microbiota Abundance and CVD Association",
    x = "Log₂ Fold Change",
    y = "–Log10 Adjusted P-Value",
    color = "Log2 FC",
    shape = "Associated with CVD"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right"
  )


########## Helper functions ###########
taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}
