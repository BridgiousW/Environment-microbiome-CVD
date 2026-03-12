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

### Panels A and B

ruralguys_withdietcategories<-read.csv('~/Novogene/codeNaturecomms/ruralguys_withdietcategories.csv')
urbanfromoriginal_subset <- read.csv('~/Novogene/urbanfromoriginal_subset.csv')

with_diet_laviiswa <- read.csv("~/Novogene/pet-table/allcombineddata.csv")

subset_data <- with_diet_laviiswa %>% select(bmi,slabno, setting, dietf,dietv,dietfh,dietm)
subset_data <- subset_data[subset_data$setting=='Rural', ]
subset_data$slabno <-paste("LS", subset_data$slabno, sep = "")

# Otu-table, CVD data, metadata - age, sex, diet  
RuralphyseqGenus <- readRDS('~/Novogene/codeNaturecomms/Ruralphyseq.rds')
#RuralphyseqGenus<-taxa_level(Ruralphyseq, "Genus")
microbiome <- subset_samples(RuralphyseqGenus, slabno%in%ruralguys_withdietcategories$slabno)
microbiome_table <- data.frame(otu_table(microbiome))

subset_metadata <- data.frame(sample_data(RuralphyseqGenus))

subset_data <- subset_data[subset_data$slabno%in%subset_metadata$slabno, ]

# Define criteria for each diet category
meat_threshold = 3        # If they eat meat on 5 or more days
fish_threshold = 5        # If they eat fish on 4 or more days
vegetarian_threshold = 3  # If they eat vegetables on 5 or more days and no meat
fruit_threshold = 3       # If they eat fruits on 5 or more days

# Categorize participants based on the criteria

subset_data <- na.omit(subset_data)
subset_data$category <- apply(subset_data[, 4:7], 1, function(row) {
  if (row["dietm"] >= meat_threshold) {
    return("Meat Lover")
  } else if (row["dietfh"] >= fish_threshold) {
    return("Fish Lover")
  } else if (row["dietv"] >= vegetarian_threshold && row["dietm"] == 0) {
    return("Vegetarian")
  } else if (row["dietf"] >= fruit_threshold) {
    return("Vegetarian")
  } else {
    return("Mixed Diet")
  }
})

ruralguys_withdietcategories$bmi <- subset_data$bmi[match(ruralguys_withdietcategories$slabno,
                                                          subset_data$slabno)]
ruralguys_withdietcategories$category <- subset_data$category[match(ruralguys_withdietcategories$slabno,
                                                          subset_data$slabno)]

ruralguys_withdietcategories <- ruralguys_withdietcategories[ruralguys_withdietcategories$slabno%in%subset_data$slabno, ]

rownames(ruralguys_withdietcategories)<-ruralguys_withdietcategories$slabno

cvddata <- ruralguys_withdietcategories[c('insulin', 'BP_Dia', 'Chol_T', 'Chol_L', 'glucose','BP_Sys','Hel_inf_status')]
confounders <- ruralguys_withdietcategories[c('age','sex','category','bmi')]

#check for relationship between CVD risk and DIET
no_of_vars = dim(cvddata)[2]
for (i in 1:no_of_vars){
  cvdrisk=cvddata[,i]
  kw <- kruskal.test(cvdrisk, confounders$category)
  print(paste(colnames(cvddata)[i], " Vs Diet +++++++++++++ pvalue = ", kw$p.value))
}

#linear regression over the microbes, confounders and cvd risk factors
# before let us make sure we have all the samples in the sample order across dataframes

infected_rural_cvddata <- cvddata[cvddata$Hel_inf_status=="Positive",]
uninfected_rural_cvddata <- cvddata[cvddata$Hel_inf_status!="Positive",]

top_microbes <- names(sort(colSums(microbiome_table), decreasing = TRUE)[1:50])
microbiome_table_top <- microbiome_table[top_microbes]

microbiome_table_top$slabno <- ruralguys_withdietcategories$slabno[
  match(rownames(microbiome_table_top), ruralguys_withdietcategories$sample.id.x)
]
microbiome_table_top <- microbiome_table_top[microbiome_table_top$slabno%in%ruralguys_withdietcategories$slabno,]

rownames(microbiome_table_top)<-microbiome_table_top$slabno
microbiome_table_top$slabno<-NULL

infected_rural_microbiome_table<-microbiome_table_top[rownames(infected_rural_cvddata), ]
uninfected_rural_microbiome_table<-microbiome_table_top[rownames(uninfected_rural_cvddata), ]

infected_rural_confounders<-confounders[rownames(infected_rural_cvddata), ]
uninfected_rural_confounders<-confounders[rownames(uninfected_rural_cvddata), ]

#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

regress <- function(cvddata, microbiome_table, confounders, group){
  # now loop over microbes
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  for (microbe in microbes){
    print(microbe)
    for (cvdrisk in cvdrisks){
      model<-lm(cvddata[,cvdrisk] ~ microbiome_table[,microbe] + confounders$age + confounders$sex + confounders$bmi + confounders$category)
      summ_model <- summary(model)
      coeffis <- summ_model$coefficients
      pvalue_ <- coeffis[2,4]
      coeffi_ <- coeffis[2,1]
      tmp<-data.frame(microbe=microbe, cvdrisk=cvdrisk, pvalue=pvalue_, coeffi_=coeffi_, group=group)
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}

infected_rural_cvddata$Hel_inf_status<-NULL
regression_infected <- regress(cvddata = infected_rural_cvddata, infected_rural_microbiome_table, infected_rural_confounders, group='Infected')

uninfected_rural_cvddata$Hel_inf_status<-NULL
regression_uninfected <- regress(uninfected_rural_cvddata, uninfected_rural_microbiome_table, uninfected_rural_confounders, group='UnInfected')

regression_combined <- rbind(regression_infected, regression_uninfected)
regression_combined$Significance <- cut(regression_combined$pvalue,
                                        breaks=c(-Inf,0.001, 0.01, 0.05, 0.1, Inf),
                                        labels = c('***', '**', '*', '', ''))

regression_combined$cvdrisk<-factor(regression_combined$cvdrisk,
                                    levels = c('BP_Sys', 'BP_Dia', 'Chol_T', 'Chol_L', 'insulin','glucose'),
                                    labels = c( "Systolic BP", "Diastolic BP",   "Total Chol", "LDL Chol", "Insulin", "Glucose"))


plot_taxa_env <- function(df){
  p <-ggplot2::ggplot(aes(x=group, y=microbe, fill=Significance), data=df)
  p <- p + ggplot2::geom_tile() + scale_fill_manual(values=rev(c("#0072B2","#F0E442","#D55E00"))) #scale_fill_b(type = 'viridis')#scale_fill_gradient2(high="#2C7BB6",low="#D7191C")
  p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
  p<-p+ ggplot2::facet_grid(. ~ cvdrisk, drop=TRUE,scale="free",space="free_x")
  #p<-p+scale_fill_steps(low="#2C7BB6",high="#D7191C", breaks = c(-Inf,0.001, 0.01, 0.05, Inf), na.value = "black")
  p<-p+ ggplot2::xlab("Groups")
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
  return(p)
}
plot_taxa_env(regression_combined)

####
# Otu-table, CVD data, metadata - age, sex, diet  
Urbanphyseq <- readRDS("~/Novogene/codeNaturecomms/Urbanphyseq.rds")
UrbanphyseqGenus<-taxa_level(Urbanphyseq, "Genus")
microbiome=subset_samples(UrbanphyseqGenus, slabno%in%urbanfromoriginal_subset$slabno)
microbiome_table <- data.frame(otu_table(microbiome))
rownames(urbanfromoriginal_subset)<- urbanfromoriginal_subset$slabno
cvddata <- urbanfromoriginal_subset[c('BP_Sys', 'BP_Dia', 'Chol_T', 'Chol_L', 'insulin','glucose','Hel_inf_status')]
confounders <- urbanfromoriginal_subset[c('age','sex')]

subset_data_urban <- with_diet_laviiswa %>% select(bmi,slabno, setting)
subset_data_urban <- subset_data_urban[subset_data_urban$setting=='Urban', ]
subset_data_urban$slabno <- sprintf("%04d", subset_data_urban$slabno)
subset_data_urban$slabno <-paste("USS", subset_data_urban$slabno, sep = "")

subset_data_urban <- subset_data_urban[subset_data_urban$slabno%in%rownames(confounders),]

confounders$bmi <- subset_data_urban$bmi[match(rownames(confounders), subset_data_urban$slabno)]

top_microbes <- names(sort(colSums(microbiome_table), decreasing = TRUE)[1:50])
microbiome_table_top <- microbiome_table[top_microbes]

microbiome_table_top$slabno <- urbanfromoriginal_subset$slabno[
  match(rownames(microbiome_table_top), urbanfromoriginal_subset$sample.id)
]
rownames(microbiome_table_top)<-microbiome_table_top$slabno
microbiome_table_top$slabno<-NULL

#linear regression over the microbes, confounders and cvd risk factors
# before let us make sure we have all the samples in the sample order across dataframes

infected_Urban_cvddata <- cvddata[cvddata$Hel_inf_status=="Positive",]
infected_Urban_cvddata$Hel_inf_status<-NULL
uninfected_Urban_cvddata <- cvddata[cvddata$Hel_inf_status!="Positive",]
uninfected_Urban_cvddata$Hel_inf_status<-NULL

infected_urban_microbiome_table<-microbiome_table_top[rownames(infected_Urban_cvddata), ]
uninfected_urban_microbiome_table<-microbiome_table_top[rownames(uninfected_Urban_cvddata), ]

infected_Urban_confounders<-confounders[rownames(infected_Urban_cvddata), ]
uninfected_urban_confounders<-confounders[rownames(uninfected_Urban_cvddata), ]

#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

regress <- function(cvddata, microbiome_table, confounders, group){
  # now loop over microbes
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  for (microbe in microbes){
    print(microbe)
    for (cvdrisk in cvdrisks){
      model<-lm(cvddata[,cvdrisk] ~ microbiome_table[,microbe] + confounders$age + confounders$sex + confounders$bmi)
      summ_model <- summary(model)
      coeffis <- summ_model$coefficients
      pvalue_ <- coeffis[2,4]
      coeffi_ <- coeffis[2,1]
      tmp<-data.frame(microbe=microbe, cvdrisk=cvdrisk, pvalue=pvalue_, coeffi_=coeffi_, group=group)
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}

regression_infected <- regress(cvddata = infected_Urban_cvddata, infected_urban_microbiome_table, infected_Urban_confounders, group='Infected')
regression_uninfected <- regress(uninfected_Urban_cvddata, uninfected_urban_microbiome_table, uninfected_urban_confounders, group='UnInfected')

regression_combined <- rbind(regression_infected, regression_uninfected)
regression_combined$Significance <- cut(regression_combined$pvalue,
                                        breaks=c(-Inf,0.001, 0.01, 0.05, 0.1, Inf),
                                        labels = c('***', '**', '*', '', ''))

regression_combined$cvdrisk<-factor(regression_combined$cvdrisk,
                                    levels = c('BP_Sys', 'BP_Dia', 'Chol_T', 'Chol_L', 'insulin','glucose'),
                                    labels = c( "Systolic BP", "Diastolic BP",   "Total Chol", "LDL Chol", "Insulin", "Glucose"))


plot_taxa_env <- function(df){
  p <-ggplot2::ggplot(aes(x=group, y=microbe, fill=Significance), data=df)
  p <- p + ggplot2::geom_tile() + scale_fill_manual(values=rev(c("#0072B2","#F0E442","#D55E00","#D7191C"))) #scale_fill_b(type = 'viridis')#scale_fill_gradient2(high="#2C7BB6",low="#D7191C")
  p<-p+ ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  p<-p+ ggplot2::geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL)
  p<-p+ ggplot2::facet_grid(. ~ cvdrisk, drop=TRUE,scale="free",space="free_x")
  #p<-p+scale_fill_steps(low="#2C7BB6",high="#D7191C", breaks = c(-Inf,0.001, 0.01, 0.05, Inf), na.value = "black")
  p<-p+ ggplot2::xlab("Groups")
  p<-p+ ggplot2::theme(strip.background = element_rect(fill = "white"))
  return(p)
}
plot_taxa_env(regression_combined)

####effect sizes
regress <- function(cvddata, microbiome_table, confounders, group){
  
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  
  for (microbe in microbes){
    for (cvdrisk in cvdrisks){
      
      # build the model
      model <- lm(
        cvddata[,cvdrisk] ~ microbiome_table[,microbe] +
          confounders$age + confounders$sex + confounders$bmi
      )
      
      summ <- summary(model)
      coeffs <- summ$coefficients
      
      # unstandardised effect size
      beta      <- coeffs[2, 1]
      pvalue    <- coeffs[2, 4]
      
      # compute standardised beta:
      # scale() converts to z-scores
      std_model <- lm(
        scale(cvddata[,cvdrisk]) ~ scale(microbiome_table[,microbe]) +
          scale(confounders$age) + confounders$sex + scale(confounders$bmi)
      )
      std_beta <- summary(std_model)$coefficients[2,1]
      
      # save results
      tmp <- data.frame(
        microbe      = microbe,
        cvdrisk      = cvdrisk,
        pvalue       = pvalue,
        beta         = beta,
        beta_std     = std_beta,
        group        = group
      )
      
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}
regression_infected <- regress(infected_Urban_cvddata,
                               infected_urban_microbiome_table,
                               infected_Urban_confounders,
                               group='Infected')

regression_uninfected <- regress(uninfected_Urban_cvddata,
                                 uninfected_urban_microbiome_table,
                                 uninfected_urban_confounders,
                                 group='UnInfected')

regression_combined_urban <- rbind(regression_infected, regression_uninfected)
#####
#### effect sizes with 95% CIs
regress <- function(cvddata, microbiome_table, confounders, group){
  
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  
  for (microbe in microbes){
    for (cvdrisk in cvdrisks){
      
      # Unstandardised model
      model <- lm(
        cvddata[,cvdrisk] ~ microbiome_table[,microbe] +
          confounders$age + confounders$sex + confounders$bmi
      )
      
      summ <- summary(model)
      coeffs <- summ$coefficients
      ci     <- confint(model)
      
      beta      <- coeffs[2, 1]
      pvalue    <- coeffs[2, 4]
      beta_lci  <- ci[2, 1]
      beta_uci  <- ci[2, 2]
      
      # Standardised model
      std_model <- lm(
        scale(cvddata[,cvdrisk]) ~ scale(microbiome_table[,microbe]) +
          scale(confounders$age) + confounders$sex + scale(confounders$bmi)
      )
      
      std_summ <- summary(std_model)
      std_ci   <- confint(std_model)
      
      beta_std     <- std_summ$coefficients[2, 1]
      beta_std_lci <- std_ci[2, 1]
      beta_std_uci <- std_ci[2, 2]
      
      # Save results
      tmp <- data.frame(
        microbe        = microbe,
        cvdrisk        = cvdrisk,
        beta           = beta,
        beta_lci       = beta_lci,
        beta_uci       = beta_uci,
        beta_std       = beta_std,
        beta_std_lci   = beta_std_lci,
        beta_std_uci   = beta_std_uci,
        pvalue         = pvalue,
        group          = group
      )
      
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}




### Panel D

require(ggalluvial)
library(MedZIM)
library(dplyr)

physeqgenus <- readRDS("~/Novogene/pet-table/physeqgenusformediation")
#otu_table(physeqgenus) <- otu_table(physeqgenus) + 1
total <- median(sample_sums(physeqgenus))
#physeqgenus <- filter_taxa(physeqgenus, function(x) sum(x > 10) > (0.10 * length(x)),TRUE)
dds<-phyloseq_to_deseq2(physeqgenus, ~ case_2)
dds<-DESeq(dds, test="Wald", fitType="local")
res1<-results(dds, cooksCutoff = FALSE)
#alpha<-0.05
#sigtab1<-res1[which(res1$padj < alpha), ]
sigtab1<-sigtab1[which(abs(sigtab1$log2FoldChange)>1), ]
sigtab1<-as.data.frame(res1)
dim(sigtab1)
#create a waterfall plot for differentially abundant taxa
df<-data.frame(sample_data(physeqgenus))["case_2"]
sigtab1$Abundant_Group <- levels(as.factor(df$case_2))[as.numeric(sigtab1$log2FoldChange>0) + 1]
df_wf<-sigtab1[c("log2FoldChange","padj", "baseMean","pvalue","Abundant_Group")]
rownames(df_wf) <- paste("taxon_", rownames(df_wf), sep = "")
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-001"] <- "taxon_Lachnospiraceae_UCG.001"
#rownames(df_wf)[rownames(df_wf) == "taxon_Lachnospiraceae_UCG-004"] <- "taxon_Lachnospiraceae_UCG.004"
#rownames(df_wf)[rownames(df_wf) == "taxon_Christensenellaceae_R-7_group"] <- "taxon_Christensenellaceae_R.7_group"

rownames(df_wf)<-gsub("-","\\.",rownames(df_wf))

mediation_df <- read.csv("~/final_mediation_schiso.csv")
med_vv <- mediation_df[mediation_df$path=="Indirect",]
filtered_med_vv <- med_vv[(med_vv$CI.2.5.. < 0 & med_vv$CI.97.5.. < 0) | (med_vv$CI.2.5.. > 0 & med_vv$CI.97.5.. > 0), ]
#head(mediation_df)
result_df <- mediation_df %>%
  filter(!(path %in% c('Direct', 'Total') & (path == 'Indirect' & pval <= 0.055)))
tmp_med <- mediation_df[mediation_df$pval<0.055,]
tmp_med <- tmp_med[tmp_med$path=="Indirect",]

result <- merge(data.frame(ID = rownames(df_wf), df_wf, row.names = NULL), tmp_med, by.x = "ID", by.y = "taxon", all.x = TRUE)
result <- merge(data.frame(ID = rownames(df_wf), df_wf, row.names = NULL), filtered_med_vv, by.x = "ID", by.y = "taxon", all.x = TRUE)

result_with_indirecteffect<- result[!is.na(result$coef), ]
colnames(result_with_indirecteffect)[colnames(result_with_indirecteffect) == "ID"] <- "taxon"
colnames(result_with_indirecteffect)[colnames(result_with_indirecteffect) == "Abundant_Group"] <- "case_2"
colnames(result_with_indirecteffect)[colnames(result_with_indirecteffect) == "case_2"] <- "Infection_status"

phylum_colors <- c(
  "red", "blue", "green", "purple", "orange",
  "cyan", "magenta", "yellow", "darkred", "darkblue",
  "darkgreen", "darkorange", "darkcyan", "darkmagenta",
  "palegoldenrod", "brown", "pink", "violet", "olive",
  "turquoise", "maroon", "navy", "gold",
  "lime", "indigo", "salmon", "peru", "orchid",
  "slateblue", "thistle", "seagreen", "darkslategray", "darkolivegreen",
  "chocolate", "firebrick", "royalblue", "darkseagreen", "darkturquoise",
  "forestgreen", "darkorchid", "indianred", "dodgerblue", "mediumvioletred",
  "orangered", "deeppink", "yellowgreen", "limegreen"
)


taxon_list=unique(tmp_med$taxon)
cvd_list<-unique(tmp_med$CVD)
tmp_med_fortable<-mediation_df[mediation_df$taxon%in%taxon_list & mediation_df$path%in%c("Indirect","Direct","Total"),]
tmp_med_fortable<-tmp_med_fortable[tmp_med_fortable$CVD%in%cvd_list,]
tmp_med_fortable[tmp_med_fortable$pval]

#tmp_med_fortable <- tmp_med[tmp_med$path%in%c("Indirect","Direct","Total"),]
tmp_med <- tmp_med[c("taxon","coef","CVD", "pval")]

ggplot(result_with_indirecteffect,
       aes(y = coef, axis1= Infection_status, axis2 = taxon, axis3 = CVD)) +
  geom_alluvium(aes(fill = taxon), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  #scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("")

###
mediation_df <- read.csv("~/final_mediation_schiso.csv")
head(mediation_df)

tmp_med <- mediation_df[mediation_df$pval<0.055,]
tmp_med <- tmp_med[tmp_med$path=="Indirect",]
tmp_med <- tmp_med[c("taxon","coef","CVD","pval")]

###
pos_med <- result_with_indirecteffect[result_with_indirecteffect$coef>0 & abs(result_with_indirecteffect$log2FoldChange)>0.2,]
pos_med$taxon<-gsub("taxon_", "", pos_med$taxon)
p<-ggplot(pos_med,aes(y = coef, axis1= Infection_status,axis2=taxon, axis3 = CVD)) +geom_alluvium(aes(fill = taxon), width = 1/12) +geom_stratum(width = 1/12,fill = "grey70", color = "black")+
  #geom_text(stat = "stratum", aes(label = after_stat(stratum), angle = 0, hjust = 0.5, vjust = 1)) +  # Adjust label angle, horizontal and vertical alignment
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_fill_manual(values = phylum_colors)+
  ggtitle("")+ theme_minimal()+  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = -.3)+ theme_minimal()+  ggrepel::geom_text_repel(
      aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA)),
      stat = "stratum", size = 5, direction = "y", nudge_x = +.3
    )+theme_void()+theme(panel.background = element_blank(),panel.border = element_rect(fill=NA),legend.text = element_text(size=12),legend.title = element_text(size = 14))
p
#ggsave("alluvial-positve_shisto_pos_mediation.pdf", p, device = "pdf", height = 7, width = 12)

#####
#####
neg_med <- result_with_indirecteffect[result_with_indirecteffect$coef<0 & abs(result_with_indirecteffect$log2FoldChange)>0.2,]
#pos_med <- result_with_indirecteffect[result_with_indirecteffect$coef>0 & abs(result_with_indirecteffect$log2FoldChange)>2.,]
neg_med$taxon<-gsub("taxon_", "", neg_med$taxon)
p<-ggplot(neg_med,
          aes(y = coef, axis1= Infection_status,axis2=taxon, axis3 = CVD)) +
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

###########get effect sizes for panel A
regress <- function(cvddata, microbiome_table, confounders, group){
  
  microbes <- names(microbiome_table)
  cvdrisks <- colnames(cvddata)
  
  df_lm_final <- data.frame()
  
  for (microbe in microbes){
    for (cvdrisk in cvdrisks){
      
      model <- lm(cvddata[,cvdrisk] ~ microbiome_table[,microbe] +
                    confounders$age + confounders$sex +
                    confounders$bmi + confounders$category)
      
      summ_model <- summary(model)
      coeffis <- summ_model$coefficients
      
      # extract p and unstandardised effect
      pvalue_ <- coeffis[2,4]
      coeffi_ <- coeffis[2,1]
      
      # compute standardised effect size
      # scale outcome + predictor
      model_std <- lm(
        scale(cvddata[,cvdrisk]) ~ scale(microbiome_table[,microbe]) +
          scale(confounders$age) + confounders$sex +
          scale(confounders$bmi) + confounders$category
      )
      std_coeff <- summary(model_std)$coefficients[2,1]
      
      # store results
      tmp <- data.frame(
        microbe = microbe,
        cvdrisk = cvdrisk,
        pvalue = pvalue_,
        beta = coeffi_,
        beta_std = std_coeff,
        group = group
      )
      
      df_lm_final <- rbind(df_lm_final, tmp)
    }
  }
  
  return(df_lm_final)
}

