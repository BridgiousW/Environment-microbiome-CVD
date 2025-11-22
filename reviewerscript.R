require(dagitty)
require(tidygraph)
require(ggdag)
setwd("/Users/mbmhvbw2/Novogene/codeNaturecomms/scripts/Rural_urban paper/Reviewerfolder")
#reviewer_mediation<-read_excel("Copy of Merged2.xlsx")
reviewer_mediation<-read.csv('allcombineddata.csv')
unique(reviewer_mediation$occpcat)

reviewer_mediation$occpcat2 <- reviewer_mediation$occpcat
reviewer_mediation$occpcat2[reviewer_mediation$occpcat%in%c('Unemployed','Housewife', 'Child/student',
                                                            'Other', 'Professional')]<-'Category 1'
reviewer_mediation$occpcat2[reviewer_mediation$occpcat%in%c('Agriculture, lumbering, charcoal',
                                                            'Fishing or lake related')]<-'Category 2'
reviewer_mediation$occpcat2[!reviewer_mediation$occpcat2%in%c('Category 1', 'Category 2')]<-'Category 3'

table(reviewer_mediation$occpcat2)

reviewer_mediation$exc4 <- reviewer_mediation$exc1
reviewer_mediation$exc4[reviewer_mediation$exc1%in%c('almost every day', 'every day')]<-'Very active'
reviewer_mediation$exc4[reviewer_mediation$exc1%in%c('less than once a month', 'not applicable','once a month')]<-'Not active'
reviewer_mediation$exc4[reviewer_mediation$exc1%in%c('once a week')]<-'A bit active'
table(reviewer_mediation$exc4)

#smoking, alcohol, occupation, exercise 
#age, sex, bmi, occp, exc1,smok,alcoh 
subset_data <- reviewer_mediation %>% dplyr::select(slabno, setting, age, sex, bmi, occpcat2, exc4,smok,alcoh)
dim(subset_data)
subset_data_urban <- subset_data[subset_data$setting=='Urban', ]
subset_data_urban$slabno <- sprintf("%04d", subset_data_urban$slabno)
subset_data_urban$slabno <-paste("USS", subset_data_urban$slabno, sep = "")
head(subset_data_urban)

subset_data_rural <- subset_data[subset_data$setting=='Rural', ]
subset_data_rural$slabno <-paste("LS", subset_data_rural$slabno, sep = "")
head(subset_data_rural)

subset_data <- rbind(subset_data_urban, subset_data_rural)

physeq <- readRDS("~/Novogene/codeNaturecomms/ruralurabnphyseq.rds")
physeq<-taxa_level(physeq, "Genus")
physeq

phy_sample_data <- data.frame(sample_data(physeq))
table(phy_sample_data$slabno%in%subset_data$slabno)

subset_data  <- subset_data[subset_data$slabno%in%phy_sample_data$slabno,]
dim(subset_data)
subset_data$infection <- phy_sample_data$Hel_inf_status[
  match(subset_data$slabno, phy_sample_data$slabno)
]
table(subset_data$infection)

subset_data$insulin <- phy_sample_data$insulin[
  match(subset_data$slabno, phy_sample_data$slabno)
]

subset_data$BP_Dia <- phy_sample_data$BP_Dia[
  match(subset_data$slabno, phy_sample_data$slabno)
]

subset_data$Chol_T <- phy_sample_data$Chol_T[
  match(subset_data$slabno, phy_sample_data$slabno)
]

subset_data$Chol_L <- phy_sample_data$Chol_L[
  match(subset_data$slabno, phy_sample_data$slabno)
]

subset_data$glucose <- phy_sample_data$glucose[
  match(subset_data$slabno, phy_sample_data$slabno)
]

subset_data$BP_Sys <- phy_sample_data$BP_Sys[
  match(subset_data$slabno, phy_sample_data$slabno)
]

dim(subset_data)
names(subset_data)

microbes <- data.frame(otu_table(physeq))
dim(microbes)

microbes$slabno <- phy_sample_data$slabno[
  match(rownames(microbes), phy_sample_data$sample.id)
]
rownames(microbes)<-microbes$slabno
microbes$slabno<-NULL
rownames(subset_data) <- subset_data$slabno

subset_data <- subset_data[rownames(microbes),]
rownames(subset_data)[1:10]
rownames(microbes)[1:10]

dim(subset_data)
dim(microbes)

# Set threshold
thresh = 0.005
microbes_filtered <- microbes[, colSums(microbes >0) / nrow(microbes) >= thresh]
dim(microbes_filtered)

####### Panel C and D ######

require(ggalluvial)
require(DESeq2)

physeqgenus <- readRDS("~/Novogene/pet-table/physeqgenusformediation")
otu_table(physeqgenus) <- otu_table(physeqgenus) + 1
total <- median(sample_sums(physeqgenus))
physeqgenus <- filter_taxa(physeqgenus, function(x) sum(x > 10) > (0.10 * length(x)),TRUE)
dds<-phyloseq_to_deseq2(physeqgenus, ~ setting.x)
dds<-DESeq(dds, test="Wald", fitType="local")
res1<-results(dds, cooksCutoff = FALSE)
alpha<-0.05
sigtab1<-res1[which(res1$padj < alpha), ]
#sigtab1<-sigtab1[which(abs(sigtab1$log2FoldChange)>1), ]
sigtab1<-as.data.frame(res1)
dim(sigtab1)
#create a waterfall plot for differentially abundant taxa
df<-data.frame(sample_data(physeqgenus))["setting.x"]
sigtab1$Abundant_Group <- levels(as.factor(df$setting.x))[as.numeric(sigtab1$log2FoldChange>0) + 1]
df_wf<-sigtab1[c("log2FoldChange","padj", "baseMean","pvalue","Abundant_Group")]
#rownames(df_wf) <- paste("taxon_", rownames(df_wf), sep = "")
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-001"] <- "taxon_Lachnospiraceae_UCG.001"
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-004"] <- "taxon_Lachnospiraceae_UCG.004"
#rownames(df_wf)[rownames(df_wf) == "taxon_Christensenellaceae_R-7_group"] <- "taxon_Christensenellaceae_R.7_group"

rownames(df_wf)<-gsub("-","\\.",rownames(df_wf))

#mediation_df <- read.csv("~/Novogene/codeNaturecomms/final_mediationRU.csv")

mediation_insulin <- read.csv('B_insulin_mediation.csv')
mediation_cholT <- read.csv('B_Chol_T_mediation.csv')
mediation_cholL <- read.csv('B_Chol_L_mediation.csv')
mediation_BPsys <- read.csv('B_BP_Sys_mediation.csv')
mediation_BPDia <- read.csv('B_BP_Dia_mediation.csv')
mediation_glucose <- read.csv('B_glucose_mediation.csv')

mediation_df <- rbind(mediation_insulin, mediation_cholT, mediation_cholL,
                      mediation_BPsys, mediation_BPDia, mediation_glucose)

sigtab1$Microbe <- rownames(sigtab1)

result <- merge.data.frame(sigtab1, mediation_df, by='Microbe')

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

###

result_with_indirecteffect$taxon <- result_with_indirecteffect$Microbe
result_with_indirecteffect$CVD <- result_with_indirecteffect$cvd_outcome

result_with_indirecteffect$CVD <- factor(
  result_with_indirecteffect$CVD,
  levels = c('insulin', 'glucose', 'BP_Sys', 'BP_Dia', 'Chol_L', 'Chol_T'),
  labels = c('Insulin', 'Glucose', 'Systolic BP', 
             'Diastolic BP', 'LDL-c', 'Total cholesterol')
)
pos_med <- result_with_indirecteffect[result_with_indirecteffect$coef>0 & abs(result_with_indirecteffect$log2FoldChange)>1,]
#pos_med <- pos_med[pos_med$pval<0.055,]

# Acidaminococcus: The indirect mediatory effect is not signifcant, and the 97% CI contains a zero
pos_med <- pos_med[pos_med$taxon!='Acidaminococcus',]

ppos<-ggplot(pos_med,
          aes(y = coef, axis1=Abundant_Group,axis2=taxon, axis3 = CVD)) +
  geom_alluvium(aes(fill = taxon), width = 1/12) +
  geom_stratum(width = 1/12,fill = "grey70", color = "black") +
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
ppos$labels$fill<-"Taxa"
ppos

#ggsave("alluvial_rural_urban_pos_mediation.pdf", ppos, device = "pdf", height = 7, width = 14)

#####
#####
neg_med <- result_with_indirecteffect[result_with_indirecteffect$coef<0 & abs(round(result_with_indirecteffect$log2FoldChange))>0.2,]

pneg<-ggplot(neg_med,
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
pneg$labels$fill<-"Taxa"
pneg


#########
library(dagitty)
library(ggdag)

# # ----- TOTAL EFFECT DAG -----
# total_dag <- dagify(
#   # outcome equations (arrows INTO node on left)
#   CVD ~ setting + microbiome + bmi + infection + exc + smok + alcohol,
#   # mediator and other nodes
#   microbiome ~ setting + age + sex + bmi + infection + exc + smok + alcohol,
#   bmi ~ setting + exc + age + sex,
#   infection ~ setting + age + sex,
#   exc ~ setting + age + sex,
#   smok ~ setting + age + sex,
#   alcohol ~ setting + age + sex,
#   # specify that age and sex are root confounders (no parents)
#   coordinates = list(x = c(setting = 0, microbiome = 1, bmi = 1.5, infection=1.5, CVD = 2),
#                      y = c(setting = 0, microbiome = 0, bmi = -1, infection = 1, CVD = 0))
# )
# 
# plot(total_dag)
# ## adjustment sets for total effect (recommended: confounders only)
# adjustmentSets(total_dag, exposure = "setting", outcome = "CVD", effect = "total")
# 
# # ----- DIRECT EFFECT DAG (block the microbiome pathway) -----
# direct_dag <- dagify(
#   CVD ~ setting + microbiome + bmi + infection + exc + smok + alcohol,
#   microbiome ~ setting + age + sex + bmi + infection + exc + smok + alcohol,
#   bmi ~ setting + exc + age + sex,
#   infection ~ setting + age + sex,
#   exc ~ setting + age + sex,
#   smok ~ setting + age + sex,
#   alcohol ~ setting + age + sex,
#   coordinates = list(x = c(setting = 0, microbiome = 1, bmi = 1.5, infection=1.5, CVD = 2),
#                      y = c(setting = 0, microbiome = 0, bmi = -1, infection = 1, CVD = 0))
# )
# 
# plot(direct_dag)
# ## adjustment sets for direct effect (this will require microbiome + confounders)
# adjustmentSets(direct_dag, exposure = "setting", outcome = "CVD", effect = "direct")
# ggdag_adjustment_set(direct_dag, exposure="setting", outcome="CVD", effect="direct", type='minimal')
# #########
# library(dagitty)
# library(ggdag)
# 
# total_dag <- dagify(
#   # outcome
#   CVD ~ setting + microbiome + bmi + infection + exc + smok + alcohol,
#   
#   # mediator + other nodes
#   microbiome ~ setting + age + sex + bmi + infection + exc + smok + alcohol,
#   bmi ~ setting + exc + age + sex,
#   infection ~ setting + age + sex,
#   exc ~ setting + age + sex,
#   smok ~ setting + age + sex,
#   alcohol ~ setting + age + sex
# )
# coordinates(total_dag) <- list(
#   x = c(setting = 0, microbiome = 1, bmi = 2, infection = 2, CVD = 3,
#         age = -1, sex = -1, exc = 1.5, smok = 1.5, alcohol = 1.5),
#   y = c(setting = 0, microbiome = 0, bmi = -1, infection = 1, CVD = 0,
#         age = -1, sex = 1, exc = -0.5, smok = 0.5, alcohol = 1.5)
# )
# plot(total_dag)
# adjustmentSets(total_dag,
#                exposure = "setting",
#                outcome = "CVD",
#                effect = "total")
#¢
total_dag <- dagify(
  CVD ~ setting + microbiome + bmi + infection + exc + smok + alcohol,
  microbiome ~ setting + age + sex + bmi + infection + exc + smok + alcohol,
  setting ~ age + sex,          
  bmi ~ setting + exc + age + sex,
  infection ~ setting + age + sex,
  exc ~ setting + age + sex,
  smok ~ setting + age + sex,
  alcohol ~ setting + age + sex
)
adjustmentSets(total_dag,
               exposure = "setting",
               outcome = "CVD",
               effect = "total")
ggdag_adjustment_set(total_dag, exposure="setting", outcome="CVD", effect="total", type='minimal')

