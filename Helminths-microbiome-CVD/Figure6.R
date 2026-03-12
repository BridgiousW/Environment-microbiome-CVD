
###### PANEL A #####

##########for total cholesterol#######
setwd("/Users/mbmhvbw2/metabolomics")
library(mixOmics)
merged_df_selected <- readRDS('~/metabolomics/merged_df_selected.rds')
t_New_fordialob_microbiome_cleaned <-read.csv("~/Novogene/codeNaturecomms/t_New_fordialob_microbiome_cleaned.csv", row.names = 1)
diabolmetadata<-read.csv("mergedmetametabolomicsimputed_2.csv", check.names = FALSE)
t_New_fordialob_microbiome_cleaned <-read.csv("~/Novogene/codeNaturecomms/t_New_fordialob_microbiome_cleaned.csv", row.names = 1)
rownames(diabolmetadata)<-diabolmetadata$slabno
diabolmetadata<-diabolmetadata[c("slabno","case")]

fortotalcholmeta <- read.csv("Totalcholesterol.csv",check.names=FALSE)
fortotalcholmeta$cluster<-gsub("^X", "",fortotalcholmeta$cluster)
fortotalcholmeta<-fortotalcholmeta[,-1]
rownames(fortotalcholmeta)<-fortotalcholmeta$cluster
fortotalcholmeta$cluster2<-gsub('.z', '/z', fortotalcholmeta$cluster)
fortotalcholmeta$HMBD_ID<-merged_df_selected$chemical_ID[match(fortotalcholmeta$cluster2, merged_df_selected$Compound)]
fortotalcholmeta$HMBD_ID[is.na(fortotalcholmeta$HMBD_ID)] <- fortotalcholmeta$cluster2[is.na(fortotalcholmeta$HMBD_ID)]
fortotalcholmeta$p_valueR<-round(unlist(as.numeric(fortotalcholmeta$p.value)),digits = 2)
fortotalcholmetaSig <- fortotalcholmeta[fortotalcholmeta$p_valueR<=0.01,]
fordialob_metabolomics <- read.csv("fordialob_metabolomics.csv",check.names=FALSE)
new_vvff_created<-t(fordialob_metabolomics)
new_vvff_created<-as.data.frame(new_vvff_created)
colnames(new_vvff_created)<-new_vvff_created[1,]
new_vvff_created<-new_vvff_created[-1,]
total_cholnew_vvff_createdSig <- new_vvff_created[gsub('/z', '.z', rownames(new_vvff_created))%in%fortotalcholmetaSig$cluster, ]
total_cholnew_vvff_createdSig$HMBD_ID<-merged_df_selected$chemical_ID[match(rownames(total_cholnew_vvff_createdSig), merged_df_selected$Compound)]
total_cholnew_vvff_createdSig$HMBD_ID[is.na(total_cholnew_vvff_createdSig$HMBD_ID)] <- row.names(total_cholnew_vvff_createdSig)[is.na(total_cholnew_vvff_createdSig$HMBD_ID)]
total_cholnew_vvff_createdSig <- total_cholnew_vvff_createdSig[!duplicated(total_cholnew_vvff_createdSig$HMBD_ID), ]
rownames(total_cholnew_vvff_createdSig)<-total_cholnew_vvff_createdSig$HMBD_ID
total_cholnew_vvff_createdSig$HMBD_ID<-NULL
total_cholnew_vvff_createdSig<-t(total_cholnew_vvff_createdSig)
total_cholnew_vvff_createdSig<-as.data.frame(total_cholnew_vvff_createdSig)
rownames_A <- rownames(total_cholnew_vvff_createdSig)
###now we get the microbes
fortotalchol_microbiome <- read.csv("totalcholmicrobe.csv",check.names=FALSE)
fortotalchol_microbiome$p_valueR<-round(unlist(as.numeric(fortotalchol_microbiome$p.value)),digits = 2)
fortotalchol_microbiomeSig <- fortotalchol_microbiome[fortotalchol_microbiome$p_valueR<=0.05,]
row.names(fortotalchol_microbiomeSig)<-fortotalchol_microbiomeSig$cluster
fortotalchol_microbiomeSig <- merge(t_New_fordialob_microbiome_cleaned,fortotalchol_microbiomeSig, by = 0, all.x = FALSE)
rownames(fortotalchol_microbiomeSig)<-fortotalchol_microbiomeSig$Row.names
fortotalchol_microbiomeSig$Row.names<-NULL
fortotalchol_microbiomeSig$cluster<-NULL
fortotalchol_microbiomeSig$p.value<-NULL
fortotalchol_microbiomeSig$Var.2<-NULL
fortotalchol_microbiomeSig$p_valueR<-NULL
fortotalchol_microbiomeSig$`cluster	`<-NULL
fortotalchol_microbiomeSig<-t(fortotalchol_microbiomeSig)
fortotalchol_microbiomeSig<-as.data.frame(fortotalchol_microbiomeSig)
##
rownames_B <- rownames(fortotalchol_microbiomeSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF <- total_cholnew_vvff_createdSig[common_row_names, ]
VVFF_num<-sapply(VVFF,as.numeric)
VVFF_num<-as.data.frame(VVFF_num)
rownames(VVFF_num)<-row.names(VVFF)
rownames_A <- rownames(VVFF_num)
rownames_B <- rownames(fortotalchol_microbiomeSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF_num <- VVFF_num[common_row_names, ]
fortotalchol_microbiomeSig <- fortotalchol_microbiomeSig[common_row_names, ]
VVFF1<-diabolmetadata[rownames(fortotalchol_microbiomeSig),]
VVFF1<-na.omit(VVFF1) 
X <- list(microbiome = fortotalchol_microbiomeSig, 
          metabolomics = VVFF_num)
Y <- VVFF1$case
summary(Y)
list.keepX <- list(microbiome = c(69, 69), metabolomics = c(93,93))
MyResult.diablo <- block.splsda(X, Y,keepX=list.keepX)
plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)
cp_totalchol<-circosPlot(MyResult.diablo, cutoff=0.70, var.adj = 0.8, size.variables = 0.7, size.labels = 0.01, color.cor = c("purple", "red"))


###### PANEL D #####

######
forldl_microbiome_editedSig <- read.csv("~/Novogene/codeNaturecomms/forldl_microbiome_editedSig.csv")
forldlometa <- read.csv("ldl.csv",check.names=FALSE)
forldlometa$cluster<-gsub("^X", "",forldlometa$cluster)
forldlometa<-forldlometa[,-1]
rownames(forldlometa)<-forldlometa$cluster
#forldlometa$cluster<-NULL
##renaming compound m/z values to their HDMBs
forldlometa$cluster2<-gsub('.z', '/z', forldlometa$cluster)
forldlometa$HMBD_ID<-merged_df_selected$chemical_ID[match(forldlometa$cluster2, merged_df_selected$Compound)]
forldlometa$HMBD_ID[is.na(forldlometa$HMBD_ID)] <- forldlometa$cluster2[is.na(forldlometa$HMBD_ID)]

#result_df$pvalue<-forldlometa$p.value[match(result_df$Compound, forldlometa$cluster)]
forldlometa$p_valueR<-round(unlist(as.numeric(forldlometa$p.value)),digits = 2)
forldlometaSig <- forldlometa[forldlometa$p_valueR<=0.01,]
fordialob_metabolomics <- read.csv("fordialob_metabolomics.csv",check.names=FALSE)
new_vvff_created<-t(fordialob_metabolomics)
new_vvff_created<-as.data.frame(new_vvff_created)
colnames(new_vvff_created)<-new_vvff_created[1,]
new_vvff_created<-new_vvff_created[-1,]
new_vvff_createdSig <- new_vvff_created[gsub('/z', '.z', rownames(new_vvff_created))%in%forldlometaSig$cluster, ]
new_vvff_createdSig$HMBD_ID<-merged_df_selected$chemical_ID[match(rownames(new_vvff_createdSig), merged_df_selected$Compound)]
new_vvff_createdSig$HMBD_ID[is.na(new_vvff_createdSig$HMBD_ID)] <- row.names(new_vvff_createdSig)[is.na(new_vvff_createdSig$HMBD_ID)]
new_vvff_createdSig <- new_vvff_createdSig[!duplicated(new_vvff_createdSig$HMBD_ID), ]
rownames(new_vvff_createdSig)<-new_vvff_createdSig$HMBD_ID
new_vvff_createdSig$HMBD_ID<-NULL
new_vvff_createdSig<-t(new_vvff_createdSig)
new_vvff_createdSig<-as.data.frame(new_vvff_createdSig)
rownames_A <- rownames(new_vvff_createdSig)
forldl_microbiome_editedSig<-t(forldl_microbiome_editedSig)
forldl_microbiome_editedSig<-as.data.frame(forldl_microbiome_editedSig)
rownames_B <- rownames(forldl_microbiome_editedSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF <- new_vvff_createdSig[common_row_names, ]
VVFF_num<-sapply(VVFF,as.numeric)
VVFF_num<-as.data.frame(VVFF_num)
rownames(VVFF_num)<-row.names(VVFF)
rownames_A <- rownames(VVFF_num)
####### microbiome data for ldl 
forldl_microbiome <- read.csv("ldlcholesterolmicrob.csv",check.names=FALSE)
forldl_microbiome$p_valueR<-round(unlist(as.numeric(forldl_microbiome$p.value)),digits = 2)
forldl_microbiome_editedSig <- forldl_microbiome[forldl_microbiome$p_valueR<=0.05,]
row.names(forldl_microbiome_editedSig)<-forldl_microbiome_editedSig$cluster
forldl_microbiome_editedSig <- merge(t_New_fordialob_microbiome_cleaned,forldl_microbiome_editedSig, by = 0, all.x = FALSE)
rownames(forldl_microbiome_editedSig)<-forldl_microbiome_editedSig$Row.names
forldl_microbiome_editedSig$Row.names<-NULL
forldl_microbiome_editedSig$cluster<-NULL
forldl_microbiome_editedSig$p.value<-NULL
forldl_microbiome_editedSig$p_valueR<-NULL
forldl_microbiome_editedSig$Var.2<-NULL
final_forldl_microbiome_editedSig<-t(forldl_microbiome_editedSig)
final_forldl_microbiome_editedSig<-as.data.frame(final_forldl_microbiome_editedSig)
rownames_B <- rownames(final_forldl_microbiome_editedSig)
common_row_names <- intersect(rownames_A, rownames_B)
VVFF_num <- VVFF_num[common_row_names, ]
final_forldl_microbiome_editedSig <- final_forldl_microbiome_editedSig[common_row_names, ]
VVFF1<-diabolmetadata[rownames(final_forldl_microbiome_editedSig),]
VVFF1<-na.omit(VVFF1) 
X <- list(microbiome = final_forldl_microbiome_editedSig, 
          metabolomics = VVFF_num)
Y <- VVFF1$case
summary(Y)
list.keepX <- list(microbiome = c(55, 55), metabolomics = c(44,44))
MyResult.diablo <- block.splsda(X, Y,keepX=list.keepX)
plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)
CP_LDL<-circosPlot(MyResult.diablo, cutoff=0.65, var.adj = 0.8, size.variables = 0.7, size.labels = 0.01, color.cor = c("purple", "red"))

##### Panel 6B ###
# Convert long to wide
corr_matrix <- dcast(cp_totalchol_long_microbes_filtered, Var1 ~ Var2, value.var = "value")
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
    title = "Totalchol-associated Microbiome–Metabolome Correlation Heatmap",
    x = "Metabolites",
    y = "Microbes",
    fill = "Correlation"
  )
#### panel 6E####
# Convert long to wide
corr_matrix <- dcast(cp_LDL_long_microbes_filtered, Var1 ~ Var2, value.var = "value")
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
    title = "LDL-associated Microbiome–Metabolome Correlation Heatmap",
    x = "Metabolites",
    y = "Microbes",
    fill = "Correlation"
  )
####panel 6c####
fortotalcholmeta$cluster2<-gsub('.z', '/z', fortotalcholmeta$cluster)
fortotalcholmeta$HMBD_ID<-merged_df_selected$chemical_ID[match(fortotalcholmeta$cluster2, merged_df_selected$Compound)]
fortotalcholmeta$HMBD_ID[is.na(fortotalcholmeta$HMBD_ID)] <- fortotalcholmeta$cluster2[is.na(fortotalcholmeta$HMBD_ID)]
fortotalcholmeta$p_valueR<-round(unlist(as.numeric(fortotalcholmeta$p.value)),digits = 2)
fortotalcholmeta <- fortotalcholmeta[fortotalcholmeta$p_valueR<=0.05,]
noduplicatefortotalcholmetasig <- fortotalcholmeta[!duplicated(fortotalcholmeta$HMBD_ID), ]
rownames(noduplicatefortotalcholmetasig)<-noduplicatefortotalcholmetasig$HMBD_ID
noduplicatefortotalcholmetasigSig3<- noduplicatefortotalcholmetasig
noduplicatefortotalcholmetasigSig3 <- noduplicatefortotalcholmetasigSig3%>%mutate(HMBD_ID = paste0("hmdb:", HMBD_ID))
###
metabolites.of.interest <- noduplicatefortotalcholmetasigSig3$HMBD_ID
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
######
