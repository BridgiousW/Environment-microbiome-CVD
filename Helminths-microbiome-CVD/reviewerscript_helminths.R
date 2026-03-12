setwd("/Users/mbmhvbw2/Novogene/codeNaturecomms/scripts/helminths paper/Reviewers")

require(ggalluvial)
require(DESeq2)
require(phyloseq)

physeqgenus <- readRDS("~/Novogene/pet-table/physeqgenusformediation")
#otu_table(physeqgenus) <- otu_table(physeqgenus) + 1
total <- median(sample_sums(physeqgenus))
#physeqgenus <- filter_taxa(physeqgenus, function(x) sum(x > 10) > (0.10 * length(x)),TRUE)
dds<-phyloseq_to_deseq2(physeqgenus, ~ case_2)
dds<-DESeq(dds, test="Wald", fitType="local")
res1<-results(dds, cooksCutoff = FALSE)
#alpha<-0.05
#sigtab1<-res1[which(res1$padj < alpha), ]
sigtab1<-as.data.frame(res1)
#sigtab1<-sigtab1[which(abs(sigtab1$log2FoldChange)>1), ]
dim(sigtab1)

rownames(sigtab1)[rownames(sigtab1) == "Lachnospiraceae_UCG-001"] <- "Lachnospiraceae_UCG.001"
rownames(sigtab1)[rownames(sigtab1) == "Lachnospiraceae_UCG-004"] <- "Lachnospiraceae_UCG.004"
rownames(sigtab1)[rownames(sigtab1) == "Christensenellaceae_R-7_group"] <- "Christensenellaceae_R.7_group"

sigtab1$taxa <- rownames(sigtab1)

#create a waterfall plot for differentially abundant taxa
df<-data.frame(sample_data(physeqgenus))["case_2"]
sigtab1$Abundant_Group <- levels(as.factor(df$case_2))[as.numeric(sigtab1$log2FoldChange>0) + 1]

# Lachnospiraceae_UCG-001, Roseburia, Christensenellaceae_R-7_group

mediation_insulin <- read.csv('insulin_helminth_mediation.csv')
mediation_cholT <- read.csv('Chol_T_helminth_mediation.csv')
mediation_cholL <- read.csv('Chol_L_helminth_mediation.csv')
mediation_BPsys <- read.csv('BP_Sys_helminth_mediation.csv')
mediation_glucose <- read.csv('glucose_helminth_mediation.csv')

mediation_df <- rbind(mediation_insulin, mediation_cholT, mediation_cholL,
                      mediation_BPsys, mediation_glucose)

sigtab1$Microbe <- rownames(sigtab1)

result <- merge.data.frame(sigtab1, mediation_df, by='Microbe')

#result <- result[result$pvalue <= 0.05,]
#write.csv(result, 'mediation-helminths.csv', quote = FALSE, row.names = TRUE)

table(result$Abundant_Group)

result$Abundant_Group <- factor(
  result$Abundant_Group,
  levels = c('Negative', 'Positive'),
  labels = c('Infected', 'Uninfected')
)

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
  levels = c('insulin', 'glucose', 'BP_Sys', 'Chol_L', 'Chol_T'),
  labels = c('Insulin', 'Glucose', 'Systolic BP', 
              'LDL-c', 'Total cholesterol')
)

pos_med <- result_with_indirecteffect[result_with_indirecteffect$coef>0 & abs(result_with_indirecteffect$log2FoldChange)>1,]

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

