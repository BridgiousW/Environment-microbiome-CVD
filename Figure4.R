
####### Panel C and D ######

require(ggalluvial)
library(MedZIM)
library(dplyr)

physeqgenus <- readRDS("~/Novogene/pet-table/physeqgenusformediation")
otu_table(physeqgenus) <- otu_table(physeqgenus) + 1
total <- median(sample_sums(physeqgenus))
physeqgenus <- filter_taxa(physeqgenus, function(x) sum(x > 10) > (0.10 * length(x)),TRUE)
dds<-phyloseq_to_deseq2(physeqgenus, ~ setting.x)
dds<-DESeq(dds, test="Wald", fitType="local")
res1<-results(dds, cooksCutoff = FALSE)
#alpha<-0.05
#sigtab1<-res1[which(res1$padj < alpha), ]
#sigtab1<-sigtab1[which(abs(sigtab1$log2FoldChange)>1), ]
sigtab1<-as.data.frame(res1)
dim(sigtab1)
#create a waterfall plot for differentially abundant taxa
df<-data.frame(sample_data(physeqgenus))["setting.x"]
sigtab1$Abundant_Group <- levels(as.factor(df$setting.x))[as.numeric(sigtab1$log2FoldChange>0) + 1]
df_wf<-sigtab1[c("log2FoldChange","padj", "baseMean","pvalue","Abundant_Group")]
rownames(df_wf) <- paste("taxon_", rownames(df_wf), sep = "")
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-001"] <- "taxon_Lachnospiraceae_UCG.001"
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-004"] <- "taxon_Lachnospiraceae_UCG.004"
#rownames(df_wf)[rownames(df_wf) == "taxon_Christensenellaceae_R-7_group"] <- "taxon_Christensenellaceae_R.7_group"

rownames(df_wf)<-gsub("-","\\.",rownames(df_wf))

mediation_df <- read.csv("~/Novogene/codeNaturecomms/final_mediationRU.csv")
head(mediation_df)
library(dplyr)

# add a taxon col on the mediation df
mediation_df <- mediation_df %>%
  mutate(
    taxon = ifelse(grepl("^taxon_", path),
                   sub("taxon_(.*?) ~ X", "\\1", path),
                   NA)
  )
# Fill down the taxon values for the remaining 4 rows in each block
mediation_df <- mediation_df %>%
  tidyr::fill(taxon, .direction = "down")
head(mediation_df)

#mediation_df <- read.csv("~/final_mediation.csv")

med_vv <- mediation_df[mediation_df$path=="Indirect",]

filtered_med_vv <- med_vv #[(med_vv$CI.2.5.. < 0 & med_vv$CI.97.5.. < 0) | (med_vv$CI.2.5.. > 0 & med_vv$CI.97.5.. > 0), ]
#head(mediation_df)
#result_df <- mediation_df %>%
#  filter(!(path %in% c('Direct', 'Total') & (path == 'Indirect' & pval <= 0.055)))
#tmp_med <- mediation_df[mediation_df$pval<0.055,]
#tmp_med <- tmp_med[tmp_med$path=="Indirect",]

#result <- merge(data.frame(ID = rownames(df_wf), df_wf, row.names = NULL), tmp_med, by.x = "ID", by.y = "taxon", all.x = TRUE)
filtered_med_vv <- filtered_med_vv[filtered_med_vv$pval<0.055,]
names(filtered_med_vv)
sigtab1$taxon <- rownames(sigtab1)

result <- merge.data.frame(sigtab1, filtered_med_vv, by='taxon')

#result <- merge(data.frame(ID = rownames(df_wf), df_wf, row.names = NULL), filtered_med_vv, by.x = "ID", by.y = "taxon", all.x = TRUE)

result_with_indirecteffect <- result
phylum_colors <- c(
  "red", "blue", "green", "purple", "orange",
  "cyan", "magenta", "yellow", "darkred", "darkblue",
  "darkgreen", "darkorange", "darkcyan", "darkmagenta",
  "palegoldenrod", "brown", "pink", "violet",
  "turquoise", "maroon", "navy", "gold",
   "indigo", "salmon", "peru", "orchid",
  "slateblue", "thistle", "seagreen", "darkslategray", "darkolivegreen",
  "chocolate", "firebrick", "royalblue", "darkseagreen", "darkturquoise",
  "forestgreen", "darkorchid", "indianred", "dodgerblue", "mediumvioletred",
  "orangered", "deeppink", "yellowgreen", "limegreen"
)


ggplot(result_with_indirecteffect,
       aes(y = coef, axis1= Abundant_Group, axis2 = taxon, axis3 = CVD)) +
  geom_alluvium(aes(fill = taxon), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("")

###
pos_med <- result_with_indirecteffect[result_with_indirecteffect$coef>0 & abs(result_with_indirecteffect$log2FoldChange)>1,]
pos_med$taxon<-gsub("taxon_", "", pos_med$taxon)
pos_med <- pos_med[pos_med$pval<0.055,]

# Acidaminococcus: The indirect mediatory effect is not signifcant, and the 97% CI contains a zero
pos_med <- pos_med[pos_med$taxon!='Acidaminococcus',]

p<-ggplot(pos_med,
          aes(y = coef, axis1=Abundant_Group,axis2=taxon, axis3 = CVD)) +
  geom_alluvium(aes(fill = taxon), width = 1/12) +
  geom_stratum(width = 1/12,fill = "grey70", color = "black") +
 # scale_fill_manual(values = phylum_colors)+
  ggtitle("")+ theme_minimal()+  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = -.3
  )+ theme_minimal()+  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = +.3
  )+theme_void()+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14))
p$labels$fill<-"Taxa"
p


#ggsave("alluvial-positve_shisto_pos_mediation.pdf", p, device = "pdf", height = 7, width = 12)

#####
#####
neg_med <- result_with_indirecteffect[result_with_indirecteffect$coef<0 & abs(round(result_with_indirecteffect$log2FoldChange))>0.2,]
neg_med$taxon<-gsub("taxon_", "", neg_med$taxon)
p<-ggplot(neg_med,
          aes(y = coef, axis1=Abundant_Group,axis2=taxon, axis3 = CVD)) +
  geom_alluvium(aes(fill = taxon), width = 1/12) +
  geom_stratum(width = 1/12,fill = "grey70", color = "black") +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum), angle = 0, hjust = 0.5, vjust = 1)) +  # Adjust label angle, horizontal and vertical alignment
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = phylum_colors)+
  ggtitle("")+ theme_minimal()+  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = -.3
  )+ theme_minimal()+  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = +.3
  )+theme_void()+
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14))
p$labels$fill<-"Taxa"
p

######### Panel A #######

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
library(reshape2)
library(ggpubr)
require(SummarizedExperiment)
require(GenomicRanges)
require(lefser)
##
ruralguys_withdietcategories<-read.csv('~/Novogene/codeNaturecomms/ruralguys_withdietcategories.csv')
urbanfromoriginal_subset <- read.csv('~/Novogene/urbanfromoriginal_subset.csv')
with_diet_laviiswa <- read.csv("~/Novogene/pet-table/allcombineddata.csv")
#subset_data <- with_diet_laviiswa %>% select(bmi,slabno, setting, dietf,dietv,dietfh,dietm)
subset_data <- dplyr::select(with_diet_laviiswa, bmi, slabno, setting, dietf, dietv, dietfh, dietm)
subset_data <- subset_data[subset_data$setting=='Rural', ]
subset_data$slabno <-paste("LS", subset_data$slabno, sep = "")
ruralguys_withdietcategories$bmi <- subset_data$bmi[match(ruralguys_withdietcategories$slabno,
                                                          subset_data$slabno)]
cvddata <- ruralguys_withdietcategories[c('slabno', 'insulin', 'BP_Dia', 'Chol_T', 'Chol_L', 'glucose','BP_Sys','Hel_inf_status')]
confounders <- ruralguys_withdietcategories[c('slabno', 'age','sex','Hel_inf_status','bmi')]
confounders$setting<-'Rural'
#Urbanphyseq <- readRDS("~/Novogene/codeNaturecomms/Urbanphyseq.rds")
#UrbanphyseqGenus<-taxa_level(Urbanphyseq, "Genus")
#microbiome=subset_samples(UrbanphyseqGenus, slabno%in%urbanfromoriginal_subset$slabno)
#microbiome_table <- data.frame(otu_table(microbiome))
rownames(urbanfromoriginal_subset)<- urbanfromoriginal_subset$slabno
cvddata1 <- urbanfromoriginal_subset[c('slabno', 'insulin', 'BP_Dia', 'Chol_T', 'Chol_L', 'glucose','BP_Sys','Hel_inf_status')]
confounders1 <- urbanfromoriginal_subset[c('slabno', 'age','sex','Hel_inf_status')]
confounders1$setting<-'Urban'

#subset_data_urban <- with_diet_laviiswa %>% select(bmi,slabno, setting)
subset_data_urban <- with_diet_laviiswa %>% dplyr::select(bmi, slabno, setting)
subset_data_urban <- subset_data_urban[subset_data_urban$setting=='Urban', ]
subset_data_urban$slabno <- sprintf("%04d", subset_data_urban$slabno)
subset_data_urban$slabno <-paste("USS", subset_data_urban$slabno, sep = "")
subset_data_urban <- subset_data_urban[subset_data_urban$slabno%in%rownames(confounders1),]
confounders1$bmi <- subset_data_urban$bmi[match(rownames(confounders1), subset_data_urban$slabno)]
# Otu-table, CVD data, metadata - age, sex, diet  
combinedforMicrobiomeCVDrisk <- rbind(cvddata, cvddata1)

physeq <- readRDS("~/Novogene/codeNaturecomms/ruralurabnphyseq.rds")
physeq<-taxa_level(physeq, "Genus")
microbiome=subset_samples(physeq, slabno%in%combinedforMicrobiomeCVDrisk$slabno)
microbiome_table <- data.frame(otu_table(microbiome))
#rownames(combinedforMicrobiomeCVDrisk)<- combinedforMicrobiomeCVDrisk$sample.id
rownames(combinedforMicrobiomeCVDrisk)<- combinedforMicrobiomeCVDrisk$slabno
cvddata <- combinedforMicrobiomeCVDrisk[c('insulin', 'BP_Dia', 'Chol_T', 'Chol_L', 'glucose', 'BP_Sys')]
confounders <- rbind(confounders, confounders1)  #combinedforMicrobiomeCVDrisk[c('age','sex','bmi','infection_status')]

meta_data <- read.delim("~/Novogene/With_case_final_novogene_received_with_lab_ids_and_outcomedata.txt", header=TRUE)
microbiome_table$slabno <- meta_data[rownames(microbiome_table), "slabno"]
rownames(microbiome_table)<-microbiome_table$slabno
microbiome_table$slabno<- NULL
##
regress <- function(cvddata, microbiome_table, confounders, group){
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  
  for (microbe in microbes){
    print(microbe)
    for (cvdrisk in cvdrisks){
      model <- lm(cvddata[[cvdrisk]] ~ microbiome_table[[microbe]] * confounders$setting +
                    confounders$age + confounders$sex + confounders$bmi + confounders$Hel_inf_status)
      
      summ_model <- summary(model)
      coeffis <- summ_model$coefficients
      
      # Extract main effect of microbe
      if ("microbiome_table[[microbe]]" %in% rownames(coeffis)) {
        pvalue_ <- coeffis["microbiome_table[[microbe]]", 4]
        coeffi_ <- coeffis["microbiome_table[[microbe]]", 1]
      } else {
        pvalue_ <- NA
        coeffi_ <- NA
      }
      
      # Extract interaction term (using grep)
      interaction_row <- grep(":", rownames(coeffis), value = TRUE)
      if (length(interaction_row) > 0 && interaction_row %in% rownames(coeffis)) {
        Interactionpvalue_ <- coeffis[interaction_row, 4]
        Interactioncoeffi_ <- coeffis[interaction_row, 1]
      } else {
        Interactionpvalue_ <- NA
        Interactioncoeffi_ <- NA
      }
      
      tmp <- data.frame(
        microbe = microbe,
        cvdrisk = cvdrisk,
        pvalue = pvalue_,
        coeffi = coeffi_,
        Interactionpvalue = Interactionpvalue_,
        Interactioncoeffi = Interactioncoeffi_,
        group = group
      )
      
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}
####
microbiome_table <- microbiome_table[rownames(cvddata),]
regression_setting <- regress(cvddata = cvddata, microbiome_table=microbiome_table,confounders = confounders, group='setting')
significant_microbes <- subset(regression_setting, pvalue <= 0.05)
significant_microbes$Inter_effect <- ifelse(significant_microbes$Interactionpvalue <= 0.05,
                                            "affected",
                                            "not affected")

library(dplyr)
library(ggplot2)
library(viridis)
####
# Clean and filter significant microbes

significant_microbes <- read.csv('~/Novogene/codeNaturecomms/significant_microbes.csv')
dim(significant_microbes)
significant_microbes_clean <- significant_microbes %>%
  filter(!is.na(Inter_effect)) %>%
  filter(Inter_effect == "not affected")

# Reorder CVD risk factors
cvd_order <- c("BP_Sys", "BP_Dia", "Chol_T", "Chol_L", "insulin", "glucose")
significant_microbes_clean$cvdrisk <- factor(significant_microbes_clean$cvdrisk, levels = cvd_order)

# Bin p-values into categories and assign asterisks + color
heatmap_data <- significant_microbes_clean %>%
  mutate(
    significance = case_when(
      pvalue <= 0.001 ~ "***",
      pvalue <= 0.01  ~ "**",
      pvalue <= 0.05  ~ "*",
      TRUE            ~ ""
    ),
    significance_color = case_when(
      significance == "***" ~ "#8B4513",  # brown
      significance == "**"  ~ "#d7191c",  # red
      significance == "*"   ~ "#2c7bb6",  # blue
      TRUE                  ~ "white"
    )
  )

# Ensure factor order of microbes to avoid overlap
heatmap_data$microbe <- factor(heatmap_data$microbe, levels = unique(heatmap_data$microbe))

# Plot

#####Not affected####
library(dplyr)
library(ggplot2)

# Clean and filter significant microbes
significant_microbes_clean<- significant_microbes %>%
  filter(!is.na(Inter_effect)) %>%
  filter(Inter_effect == "not affected")

# Reorder CVD risk factors
#cvd_order <- c("sbpnew", "dbpnew", "cholesterol2", "ldl_c", "insulin", "glucose")
#significant_microbes_clean$cvdrisk <- factor(significant_microbes_clean$cvdrisk, levels = cvd_order)

# Bin p-values into categories for coloring and asterisks
heatmap_data <- significant_microbes_clean %>%
  mutate(
    significance = case_when(
      pvalue <= 0.001 ~ "***",
      pvalue <= 0.01  ~ "**",
      pvalue <= 0.05  ~ "*",
      TRUE            ~ ""
    ),
    pval_class = case_when(
      pvalue <= 0.001 ~ "p ≤ 0.001",
      pvalue <= 0.01  ~ "p ≤ 0.01",
      pvalue <= 0.05  ~ "p ≤ 0.05",
      TRUE            ~ NA_character_
    )
  ) %>%
  filter(!is.na(pval_class))  # Keep only significant ones

# Ensure proper factor levels
heatmap_data$microbe <- factor(heatmap_data$microbe, levels = unique(heatmap_data$microbe))
heatmap_data$pval_class <- factor(heatmap_data$pval_class, levels = c("p ≤ 0.05", "p ≤ 0.01", "p ≤ 0.001"))

# Plot
ggplot(heatmap_data, aes(x = cvdrisk, y = microbe)) +
  geom_tile(aes(fill = pval_class), color = "white") +
  geom_text(aes(label = significance), color = "black", size = 6, family = "Helvetica", fontface = "bold") +
  scale_fill_manual(
    name = "Significance level",
    values = c(
      "p ≤ 0.05" = "#2c7bb6",  # blue
      "p ≤ 0.01" = "#d7191c",  # red
      "p ≤ 0.001" = "#8B4513"  # brown
    )
  ) +
  labs(
    title = "Strength of Microbe–CVD Risk Associations (Not Affected by Setting)",
    x = "CVD Risk Factor",
    y = "Microbe"
  ) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


#### Panel B #####
#####affected by setting
library(ggplot2)
library(ggrepel)
library(dplyr)

# Clean up names
sigtab <- read.csv("~/Novogene/pet-table/sigtab.csv", row.names = 1)
sigtab$taxa <- rownames(sigtab)
sigtab$taxa_clean <- gsub("_", " ", sigtab$taxa)
significant_microbes$microbe_clean <- gsub("_", " ", significant_microbes$microbe)

# Mark overlap and enrichment group
sigtab <- sigtab %>%
  mutate(
    overlap_with_cvd = ifelse(taxa_clean %in% significant_microbes$microbe_clean, "CVD-associated", "Not associated"),
    highlight_group = case_when(
      overlap_with_cvd == "CVD-associated" & log2FoldChange > 0 ~ "CVD-associated (Urban)",
      overlap_with_cvd == "CVD-associated" & log2FoldChange < 0 ~ "CVD-associated (Rural)",
      TRUE ~ "Not associated"
    )
  )

# Filter for labeling: top 25 overlapping taxa by lowest padj
top_labels <- sigtab %>%
  filter(overlap_with_cvd == "CVD-associated") %>%
  arrange(padj) %>%
  slice(1:25)

# Volcano plot
ggplot(sigtab, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = highlight_group), alpha = 0.8, size = 2.5) +
  geom_text_repel(
    data = top_labels,
    aes(label = paste0("italic('", taxa_clean, "')")),
    parse = TRUE,
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40"
  ) +
  scale_color_manual(values = c(
    "CVD-associated (Urban)" = "#1b9e77",
    "CVD-associated (Rural)" = "#d95f02",
    "Not associated" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: Differentially Abundant Taxa Associated with CVD Risk",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Microbe Status"
  )

