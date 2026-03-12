require(mixOmics)

# Panel D
##systolicfinalwithHDMI
diabolmetadata<-read.csv("mergedmetametabolomicsimputed_2.csv", check.names = FALSE)
t_New_fordialob_microbiome_cleaned <-read.csv("~/Novogene/codeNaturecomms/t_New_fordialob_microbiome_cleaned.csv", row.names = 1)
rownames(diabolmetadata)<-diabolmetadata$slabno
diabolmetadata<-diabolmetadata[c("slabno","case")]
fordialob_metabolomics <- read.csv("fordialob_metabolomics.csv",check.names=FALSE)
systol_Sig_microbiome <- read.csv("systolicmicrob.csv",check.names=FALSE)
forsystolicmeta <- read.csv("sysBp.csv",check.names=FALSE)
forsystolicmeta$cluster<-gsub("^X", "",forsystolicmeta$cluster)
forsystolicmeta<-forsystolicmeta[,-1]
rownames(forsystolicmeta)<-forsystolicmeta$cluster
forsystolicmeta$p_valueR<-round(unlist(as.numeric(forsystolicmeta$p.value)),digits = 2)
forsystolicmetaSig <- forsystolicmeta[forsystolicmeta$p_valueR<=0.01,]
###second attempt SYS using VVCONHUGATION
forsystolicmetaSig$cluster2<-gsub('.z', '/z', forsystolicmetaSig$cluster)
forsystolicmetaSig$HMBD_ID<-merged_df_selected$chemical_ID[match(forsystolicmetaSig$cluster2, merged_df_selected$Compound)]
forsystolicmetaSig$HMBD_ID[is.na(forsystolicmetaSig$HMBD_ID)] <- forsystolicmetaSig$cluster2[is.na(forsystolicmetaSig$HMBD_ID)]
forsystolicmetaSig$p_valueR<-round(unlist(as.numeric(forsystolicmetaSig$p.value)),digits = 2)
forsystolicmetaSig <- forsystolicmetaSig[forsystolicmetaSig$p_valueR<=0.05,]
new_vvff_created<-t(fordialob_metabolomics)
new_vvff_created<-as.data.frame(new_vvff_created)
colnames(new_vvff_created)<-new_vvff_created[1,]
new_vvff_created<-new_vvff_created[-1,]
systol_new_vvff_createdSig <- new_vvff_created[gsub('/z', '.z', rownames(new_vvff_created))%in%forsystolicmetaSig$cluster, ]
systol_new_vvff_createdSig$HMBD_ID<-merged_df_selected$chemical_ID[match(rownames(systol_new_vvff_createdSig), merged_df_selected$Compound)]
systol_new_vvff_createdSig$HMBD_ID[is.na(systol_new_vvff_createdSig$HMBD_ID)] <- row.names(systol_new_vvff_createdSig)[is.na(systol_new_vvff_createdSig$HMBD_ID)]
systol_new_vvff_createdSig <- systol_new_vvff_createdSig[!duplicated(systol_new_vvff_createdSig$HMBD_ID), ]
rownames(systol_new_vvff_createdSig)<-systol_new_vvff_createdSig$HMBD_ID
systol_new_vvff_createdSig$HMBD_ID<-NULL
systol_new_vvff_createdSig<-t(systol_new_vvff_createdSig)
systol_new_vvff_createdSig<-as.data.frame(systol_new_vvff_createdSig)
rownames_A <- rownames(systol_new_vvff_createdSig)
#result_df$pvalue<-forldlometa$p.value[match(result_df$Compound, forldlometa$cluster)]
###now we get the microbes
systol_Sig_microbiome$p_valueR<-round(unlist(as.numeric(systol_Sig_microbiome$p.value)),digits = 2)
systol_Sig_microbiome <- systol_Sig_microbiome[systol_Sig_microbiome$p_valueR<=0.05,]
row.names(systol_Sig_microbiome)<-systol_Sig_microbiome$cluster
systol_Sig_microbiome <- merge(t_New_fordialob_microbiome_cleaned,systol_Sig_microbiome, by = 0, all.x = FALSE)
rownames(systol_Sig_microbiome)<-systol_Sig_microbiome$Row.names
systol_Sig_microbiome$Row.names<-NULL
systol_Sig_microbiome$cluster<-NULL
systol_Sig_microbiome$p.value<-NULL
systol_Sig_microbiome$Var.2<-NULL
systol_Sig_microbiome$p_valueR<-NULL
systol_Sig_microbiome$`cluster	`<-NULL
systol_Sig_microbiome<-t(systol_Sig_microbiome)
systol_Sig_microbiome<-as.data.frame(systol_Sig_microbiome)
##
rownames_B <- rownames(systol_Sig_microbiome)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF <- systol_new_vvff_createdSig[common_row_names, ]
VVFF_num<-sapply(VVFF,as.numeric)
VVFF_num<-as.data.frame(VVFF_num)
rownames(VVFF_num)<-row.names(VVFF)
rownames_A <- rownames(VVFF_num)
rownames_B <- rownames(systol_Sig_microbiome)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF_num <- VVFF_num[common_row_names, ]
systol_Sig_microbiome <- systol_Sig_microbiome[common_row_names, ]
VVFF1<-diabolmetadata[rownames(systol_Sig_microbiome),]
VVFF1<-na.omit(VVFF1) 
X <- list(microbiome = systol_Sig_microbiome, 
          metabolomics = VVFF_num)
Y <- VVFF1$case
summary(Y)
list.keepX <- list(microbiome = c(63, 63), metabolomics = c(42,42))
MyResult.diablo <- block.splsda(X, Y,keepX=list.keepX)
plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)
cp_sys<-circosPlot(MyResult.diablo, cutoff=0.70, var.adj = 0.8, size.variables = 0.7, size.labels = 0.01, color.cor = c("purple", "red"))

# Panel A
######for diastolic with HDMI
fordiastolicmeta <- read.csv("diabp.csv",check.names=FALSE)
fordiastolicmeta$cluster<-gsub("^X", "",fordiastolicmeta$cluster)
fordiastolicmeta<-fordiastolicmeta[,-1]
rownames(fordiastolicmeta)<-fordiastolicmeta$cluster
fordiastolicmeta$cluster2<-gsub('.z', '/z', fordiastolicmeta$cluster)
fordiastolicmeta$HMBD_ID<-merged_df_selected$chemical_ID[match(fordiastolicmeta$cluster2, merged_df_selected$Compound)]
fordiastolicmeta$HMBD_ID[is.na(fordiastolicmeta$HMBD_ID)] <- fordiastolicmeta$cluster2[is.na(fordiastolicmeta$HMBD_ID)]
fordiastolicmeta$p_valueR<-round(unlist(as.numeric(fordiastolicmeta$p.value)),digits = 2)
fordiastolicmetaSig <- fordiastolicmeta[fordiastolicmeta$p_valueR<=0.01,]
fordialob_metabolomics <- read.csv("fordialob_metabolomics.csv",check.names=FALSE)
new_vvff_created<-t(fordialob_metabolomics)
new_vvff_created<-as.data.frame(new_vvff_created)
colnames(new_vvff_created)<-new_vvff_created[1,]
new_vvff_created<-new_vvff_created[-1,]
fordiastolicmetaSig <- new_vvff_created[gsub('/z', '.z', rownames(new_vvff_created))%in%fordiastolicmetaSig$cluster, ]
fordiastolicmetaSig$HMBD_ID<-merged_df_selected$chemical_ID[match(rownames(fordiastolicmetaSig), merged_df_selected$Compound)]
fordiastolicmetaSig$HMBD_ID[is.na(fordiastolicmetaSig$HMBD_ID)] <- row.names(fordiastolicmetaSig)[is.na(fordiastolicmetaSig$HMBD_ID)]
fordiastolicmetaSig <- fordiastolicmetaSig[!duplicated(fordiastolicmetaSig$HMBD_ID), ]
rownames(fordiastolicmetaSig)<-fordiastolicmetaSig$HMBD_ID
fordiastolicmetaSig$HMBD_ID<-NULL
fordiastolicmetaSig<-t(fordiastolicmetaSig)
fordiastolicmetaSig<-as.data.frame(fordiastolicmetaSig)
rownames_A <- rownames(fordiastolicmetaSig)
###now we get the microbes
fordiastolic_microbiome <- read.csv("diastolicmicrob.csv",check.names=FALSE)
fordiastolic_microbiome$p_valueR<-round(unlist(as.numeric(fordiastolic_microbiome$p.value)),digits = 2)
fordiastolic_microbiomeSig <- fordiastolic_microbiome[fordiastolic_microbiome$p_valueR<=0.05,]
row.names(fordiastolic_microbiomeSig)<-fordiastolic_microbiomeSig$cluster
fordiastolic_microbiomeSig <- merge(t_New_fordialob_microbiome_cleaned,fordiastolic_microbiomeSig, by = 0, all.x = FALSE)
rownames(fordiastolic_microbiomeSig)<-fordiastolic_microbiomeSig$Row.names
fordiastolic_microbiomeSig$Row.names<-NULL
fordiastolic_microbiomeSig$cluster<-NULL
fordiastolic_microbiomeSig$p.value<-NULL
fordiastolic_microbiomeSig$Var.2<-NULL
fordiastolic_microbiomeSig$p_valueR<-NULL
fordiastolic_microbiomeSig$`cluster	`<-NULL
fordiastolic_microbiomeSig<-t(fordiastolic_microbiomeSig)
fordiastolic_microbiomeSig<-as.data.frame(fordiastolic_microbiomeSig)
##
rownames_B <- rownames(fordiastolic_microbiomeSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF <- fordiastolicmetaSig[common_row_names, ]
VVFF_num<-sapply(VVFF,as.numeric)
VVFF_num<-as.data.frame(VVFF_num)
rownames(VVFF_num)<-row.names(VVFF)
rownames_A <- rownames(VVFF_num)
rownames_B <- rownames(fordiastolic_microbiomeSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF_num <- VVFF_num[common_row_names, ]
fordiastolic_microbiomeSig <- fordiastolic_microbiomeSig[common_row_names, ]
VVFF1<-diabolmetadata[rownames(fordiastolic_microbiomeSig),]
VVFF1<-na.omit(VVFF1) 
X <- list(microbiome = fordiastolic_microbiomeSig, 
          metabolomics = VVFF_num)
Y <- VVFF1$case
summary(Y)
list.keepX <- list(microbiome = c(97,97), metabolomics = c(118,118))
MyResult.diablo <- block.splsda(X, Y,keepX=list.keepX)
plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)
cp_diast<-circosPlot(MyResult.diablo, cutoff=0.70, var.adj = 0.8, size.variables = 0.7, size.labels = 0.01, color.cor = c("purple", "red"))

#####Panel 7B
# Convert long to wide
corr_matrix <- dcast(cp_dias_long_microbes_filtered, Var1 ~ Var2, value.var = "value")
rownames(corr_matrix) <- corr_matrix$Var1
corr_matrix <- corr_matrix[, -1]

# Transpose to switch axes: Microbes on Y, Metabolites on X
corr_matrix <- t(as.matrix(corr_matrix))

# Set exact 0 to NA so it appears white
corr_matrix[corr_matrix == 0] <- NA

# Heatmap
ggcorrplot(
  corr_matrix,
  lab = FALSE,
  type = "full",
  ggtheme = theme_minimal(),
  colors = c("#2166AC", "white", "#B2182B"),
  outline.col = "white"
) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "white",
    limits = c(-1, 1)
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "Diastolic-associated Microbiome–Metabolome Correlation Heatmap",
    x = "Metabolites",
    y = "Microbes",
    fill = "Correlation"
  )
####
######## Panel 7c #####
fordiastolicmetaSig$HMBD_ID<-merged_df_selected$chemical_ID[match(fordiastolicmetaSig$cluster2, merged_df_selected$Compound)]
fordiastolicmetaSig$HMBD_ID[is.na(fordiastolicmetaSig$HMBD_ID)] <- fordiastolicmetaSig$cluster2[is.na(fordiastolicmetaSig$HMBD_ID)]
noduplicatefordiastolicmetaSig <- fordiastolicmetaSig[!duplicated(fordiastolicmetaSig$HMBD_ID), ]
rownames(noduplicatefordiastolicmetaSig)<-noduplicatefordiastolicmetaSig$HMBD_ID
noduplicatefordiastolicmetaSig3<- noduplicatefordiastolicmetaSig
noduplicatefordiastolicmetaSig3 <- noduplicatefordiastolicmetaSig3%>%mutate(HMBD_ID = paste0("hmdb:", HMBD_ID))
###
metabolites.of.interest <- noduplicatefordiastolicmetaSig3$HMBD_ID
chemical.classes <- chemicalClassSurvey(db = rampDB, mets = metabolites.of.interest)
metabolite.classes <- as.data.frame(chemical.classes$met_classes)
#datatable(metabolite.classes)
df_classes <- metabolite.classes
ClassyFire_class <- df_classes[df_classes$class_level_name=="ClassyFire_class",]
class_firedf <- data.frame(table(ClassyFire_class$class_name))
head(class_firedf)
class_firedf <- class_firedf[class_firedf$Freq>1,]
class_firedf <- class_firedf[order(class_firedf$Freq),]
class_firedf$Var1<-factor(class_firedf$Var1, levels = class_firedf$Var1)
ggplot(data = class_firedf,
       aes(x=Var1,y=Freq))+
  geom_bar(stat='identity',fill='#6082B6')+coord_flip()+ylab('No. of metabolites')+xlab('')+
  theme(axis.title = element_text(color="black", size=14))+
  theme(axis.text = element_text(color="black", size=12))+
  theme(panel.background = element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.background = element_rect(fill = "white",colour = "black"),
        strip.text = element_text(size=12),
        panel.border = element_rect(fill=NA))+ggtitle(label = 'ClassyFire_class')
####Panel 7E#####
# Convert long to wide
corr_matrix <- dcast(cp_cp_sys_long_microbes_filtered, Var1 ~ Var2, value.var = "value")
rownames(corr_matrix) <- corr_matrix$Var1
corr_matrix <- corr_matrix[, -1]

# Transpose to switch axes: Microbes on Y, Metabolites on X
corr_matrix <- t(as.matrix(corr_matrix))

# Set exact 0 to NA so it appears white
corr_matrix[corr_matrix == 0] <- NA

# Heatmap
ggcorrplot(
  corr_matrix,
  lab = FALSE,
  type = "full",
  ggtheme = theme_minimal(),
  colors = c("#2166AC", "white", "#B2182B"),
  outline.col = "white"
) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "white",
    limits = c(-1, 1)
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "systolic BP-associated Microbiome–Metabolome Correlation Heatmap",
    x = "Metabolites",
    y = "Microbes",
    fill = "Correlation"
  )
