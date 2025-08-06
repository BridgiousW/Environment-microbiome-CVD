###Panel A #####

# Helper functions
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

#####Diffeabundanceruralurban.
tax_lev<-"Genus"
physeq<-readRDS("~/Novogene/codeNaturecomms/PhyloseqObjectHelminths.rds")
physeqcasegenus<-taxa_level(physeq, tax_lev)
otu_table(physeqcasegenus)<-otu_table(physeqcasegenus)+1
dds<-phyloseq_to_deseq2(physeqcasegenus, ~ setting.x)
#dds <- phyloseq_to_deseq2(physeqcasegenus, ~ age + sex + bmi + Hel_inf_status + setting.x)
dds<-DESeq(dds, test="Wald", fitType="local")

res1<-results(dds, cooksCutoff = FALSE)
alpha<-0.05
sigtab<-res1[which(res1$padj < alpha), ]
sigtab<-sigtab[which(abs(sigtab$log2FoldChange)>=1.0), ]
sigtab<-as.data.frame(sigtab)
#create a waterfall plot for differentially abundant taxa
df<-data.frame(sample_data(physeqcasegenus))["setting.x"]
sigtab$Abundant_Group <- levels(as.factor(df$setting.x))[as.numeric(sigtab$log2FoldChange>0) + 1]
####
library(dplyr)
library(ggplot2)

# Filter top 80 taxa by absolute log2FoldChange
top_taxa <- sigtab %>%
  arrange(-abs(log2FoldChange)) %>%
  head(80)

# Clean and reorder taxa names
top_taxa$taxa <- rownames(top_taxa)
top_taxa$taxa <- gsub("_", " ", top_taxa$taxa)
top_taxa$taxa <- factor(top_taxa$taxa, levels = rev(top_taxa$taxa))

# Plot
ggplot(top_taxa, aes(x = taxa, y = log2FoldChange, fill = Abundant_Group)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("Rural" = "#1b9e77", "Urban" = "#d95f02")) +
  labs(
    title = "Top 80 Differentially Abundant Genera",
    x = NULL, y = "Log2 Fold Change",
    fill = "Enriched In"
  ) +
  theme(
    axis.text.y = element_text(size = 7, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey85")
  )

###### ###Panel B #####
res_unique <- read.csv('~/Novogene/pet-table/res_unique.csv')
#sigtab$Genus <- rownames(sigtab)
sigtab2<-read.csv("~/Novogene/pet-table/sigtab.csv")
# Assign 'Urban' if score > 0, else 'Rural'
res_unique$setting_enriched <- ifelse(res_unique$scores > 0, "Urban", "Rural")

#####
#install.packages("VennDiagram")
library(VennDiagram)
# Get sets of genera enriched in each setting by each method

# DESeq2 results
#urban_deseq  <- sigtab$Genus[sigtab$Abundant_Group == "Urban"]
#rural_deseq  <- sigtab$Genus[sigtab$Abundant_Group == "Rural"]

urban_deseq  <- sigtab2$Genus[sigtab2$Abundant_Group == "Urban"]
rural_deseq  <- sigtab2$Genus[sigtab2$Abundant_Group == "Rural"]

# LEfSe results
urban_lefse  <- res_unique$Genus[res_unique$setting_enriched == "Urban"]
rural_lefse  <- res_unique$Genus[res_unique$setting_enriched == "Rural"]

###
venn.plot <- venn.diagram(
  x = list(
    "Urban (DESeq2)" = unique(urban_deseq),
    "Rural (DESeq2)" = unique(rural_deseq),
    "Urban (LEfSe)"  = unique(urban_lefse),
    "Rural (LEfSe)"  = unique(rural_lefse)
  ),
  filename = NULL,  # NULL so we can plot to screen
  fill = c("red", "green", "blue", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.col = "black",
  margin = 0.05,
  main = "Overlap of Enriched Genera from LDA and DA"
)

grid::grid.draw(venn.plot)
