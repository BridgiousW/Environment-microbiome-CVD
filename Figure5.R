#Load the nlme package (and other relevant packages)
library(ggplot2)
library(nlme)
library(lmtest)
library(dplyr)
library(knitr)
library(readxl)
library(writexl)
library(simpleFDR)
library(fuzzySim)
library(ggrepel)
library(ggplot2)
library(dplyr)
####Panel A######
setwd("/Users/mbmhvbw2/metabolomics")
#Load meta file and create data frame(s)
meta <- read.csv("base2transwithout194_overall.csv") #load the meta file
#####
meta$setting <- factor(meta$setting, levels = c("Rural", "Urban"))
(metaB <- as.data.frame(meta[, c(52:6559)]))
df = data.frame()
for (i in 1:ncol(metaB)) { 
  meta_name <- names(metaB)[i]
  model <- lm(metaB[, i] ~ setting + age.x + sex.x, data = meta)
  coef_summary <- summary(model)$coefficients
  
  # Use grep to locate the row for 'urban'
  row_name <- grep("Urban", rownames(coef_summary), value = TRUE)
  
  if (length(row_name) > 0) {
    log2FC <- coef_summary[row_name, 1]
    p_value <- coef_summary[row_name, 4]
    
    tmp <- data.frame(
      Metabolite = meta_name,
      log2FC = log2FC,
      p_value = p_value
    )
    df <- rbind(df, tmp)
  }
}
df$Metabolite<-gsub("^X", "",df$Metabolite)
df$Metabolite<-gsub('.z', '/z', df$Metabolite)
df$HMBD_ID<-merged_df_selected$chemical_ID[match(df$Metabolite, merged_df_selected$Compound)]
df$HMBD_ID[is.na(df$HMBD_ID)] <- df$Metabolite[is.na(df$HMBD_ID)]
# Add a significance column
df <- df %>%
  mutate(significance = case_when(
    p_value < 0.01 & log2FC > 1  ~ "Up in Urban",
    p_value < 0.01 & log2FC < -1 ~ "Down in Urban",
    TRUE                         ~ "Not Significant"
  ))
####
top_hits <- df[df$p_value < 0.01 & abs(df$log2FC) > 1, ]
top_hits <- top_hits[order(top_hits$p_value), ]
top_hits <- head(top_hits, 20)  # label top 10

ggplot(df, aes(x = log2FC, y = -log10(p_value), color = significance)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_text_repel(
    data = top_hits,
    aes(label = HMBD_ID),  # change to your metabolite name column
    size = 4,
    max.overlaps = 15,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50"
  ) +
  scale_color_manual(
    values = c("Up in Urban" = "#D7191C",     # red
               "Down in Urban" = "#2C7BB6",   # blue
               "Not Significant" = "grey70")  # grey
  ) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Metabolite Differences (Urban vs Rural)",
    x = expression(Log[2]~Fold~Change~(Urban~vs~Rural)),
    y = expression(-Log[10]~p~value),
    color = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),       # Removes all grid lines
    panel.background = element_blank(), # Removes background color
    plot.background = element_blank()   # Removes outer plot background
  )
###
#####Panel B #######

library(RaMP)
rampDB <- RaMP()
top_hitsfoeIMPALA <- df[df$p_value < 0.05 & abs(df$log2FC) > 1, ]
HMDB_top_hitsfoeIMPALA<-top_hitsfoeIMPALA$HMBD_ID
HMDB_top_hitsfoeIMPALA<-as.data.frame(HMDB_top_hitsfoeIMPALA)
hmdb_only <- HMDB_top_hitsfoeIMPALA[grepl("^HMDB\\d+$", HMDB_top_hitsfoeIMPALA$HMDB_top_hitsfoeIMPALA), ]
hmdb_only<-as.data.frame(hmdb_only)
##loop over metabolite names.
hmdb_only <- hmdb_only%>%mutate(HMBD_ID = paste0("hmdb:", HMBD_ID))
metabolites.of.interest <- hmdb_only$HMBD_ID
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
#####Panel D #####

fisher.results <- runCombinedFisherTest(db = rampDB, analytes = metabolites.of.interest)
datatable(data.frame(fisher.results))
filtered.fisher.results <- FilterFishersResults(fisher.results, pval_type = 'holm', pval_cutoff=0.05)
clusters <- RaMP::findCluster(db = rampDB, filtered.fisher.results,
                              perc_analyte_overlap = 0.2,
                              min_pathway_tocluster = 2, perc_pathway_overlap = 0.2
)
datatable(clusters$fishresults %>% mutate_if(is.numeric, ~ round(., 8)),
          rownames = FALSE
)
pathwayResultsPlot(db = rampDB, fisher.results, text_size = 10, perc_analyte_overlap = 0.2, 
                   min_pathway_tocluster = 2, perc_pathway_overlap = 0.2, interactive = FALSE)


rdmp_ldl_pathways<-getPathwayFromAnalyte(db = rampDB,metabolites.of.interest)
df <- rdmp_ldl_pathways[c('pathwaySource','pathwayName')]
df_tab <- data.frame(table(df$pathwaySource, df$pathwayName))
head(df_tab)

df_tab <- df_tab[df_tab$Freq>0,]
####Panel C#####
df<- read.csv("DataforPLSDA.csv",check.names=FALSE)
plsdamodel<-plslda(setting~., data = df, pls="oscorespls")
#summary(plsdamodel)
summary(plsdamodel$pls.out)

comp1 <- plsdamodel$pls.out$scores[,1]
comp2 <- plsdamodel$pls.out$scores[,2]

df2 <- data.frame(comp1, comp2, Setting=df$setting)

p <- ggplot(data = df2, aes(x=comp1, y=comp2, color=Setting))+
  geom_point()+theme(axis.title.x = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.x = element_text(color="black", size=12))+
  theme(axis.title.y = element_text(color="black", size=15,face = "bold"))+
  theme(axis.text.y = element_text(color="black", size=12))+
  scale_color_manual("Setting",values=rev(c("#DAA520","#4682B4"))) +
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

######Panel E ####
library(ropls)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")
library(pROC)
install.packages("ropls")
setwd("/Users/mbmhvbw2/metabolomics")
# Generate a random example dataset with 100 metabolites and 100 samples
set.seed(42)  # Set a seed for reproducibility
num_metabolites <- 6508
num_samples <- 207


DataforPLSDA<- read.csv("DataforPLSDA.csv",check.names=FALSE)



# Generate metabolite data
metabolite_data <- matrix(DataforPLSDA[,8:6515])
#metabolite_data <- DataforPLSDA[,8:6515]

############ LETS CHECK WITH TRAIN SPLIT#############

require(caret)
indexMetabols<-createDataPartition(DataforPLSDA$setting, p=0.7, list = FALSE)
trainData <- DataforPLSDA[indexMetabols, ]
testData <- DataforPLSDA[-indexMetabols, ]

table(trainData$setting)

exclude_vars <- c(1, 3:7)
exclude_vars

trainData <- trainData[-exclude_vars]

ctrl1 <- caret::trainControl(method = "repeatedcv", 
                             repeats = 3, number = 5, 
                             classProbs = TRUE,
                             savePredictions=TRUE,
                             summaryFunction=twoClassSummary)
ctrl2 <- caret::trainControl(method = "cv", 
                             number = 5, 
                             classProbs = TRUE,
                             savePredictions=TRUE,
                             summaryFunction=twoClassSummary)

model1 <- caret::train(setting~.,
                       data=trainData,
                       method="pls",
                       trControl=ctrl1)
model2 <- caret::train(setting~.,
                       data=trainData,
                       method="svm",
                       trControl=ctrl2)

model1
###

# Generate class labels
confusionMatrix(predict(model1, testData), as.factor(testData$setting))

data_pred <- as.numeric(model1$pred$pred)
data_test <- as.numeric(model1$pred$obs)
roc_ <- roc(data_test, data_pred)
roc_

ci.auc(roc_)

gbImp <- varImp(model1, scale = FALSE)
gbImp  

pls.pobs <- predict(model1, testData, type="prob")
head(pls.pobs)

###
pls.ROC <- roc(predictor=pls.pobs$Urban,
               response=testData$setting)
#levels=rev(levels(as.factor(testData$setting))))
plot(pls.ROC,main="pls-roc")
ggroc(pls.ROC) + geom_abline(intercept = 1,slope = 1)+ylim(c(0,1))+xlim(c(1,0))
####
g <- ggroc(pls.ROC)
g
# with additional aesthetics:
k<-ggroc(pls.ROC, alpha = 0.5, colour = "blue", linetype = 1, size = 1)
k
k+ theme_minimal() + ggtitle("ROC curve showing accuracy of the metabolites-based predictive model ") +  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="red", linetype="dashed")
# You can then your own theme, etc.
k + theme_minimal() + ggtitle("ROC curve showing accuracy of the metabolites-based predictive model ") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="red", , linetype="dashed")
#####
