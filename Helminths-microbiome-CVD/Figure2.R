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

#PANEL A
#import novogene alpha diversity indices
observed.diversity <- read.delim("~/Novogene/core-metrics-results/observed_features_vector.tsv/alpha-diversity.tsv", row.names = 1)
observed.diversity$sample<-rownames(observed.diversity)
shannon.diversity <- read.delim("~/Novogene/core-metrics-results/shannon_vector.tsv/alpha-diversity.tsv", row.names = 1)
shannon.diversity$sample<-rownames(shannon.diversity)
meta_data <- read.delim("~/Novogene/With_case_final_novogene_received_with_lab_ids_and_outcomedata.txt", header=TRUE)
rownames(meta_data)<-meta_data$sample.id
meta_data_2 <- meta_data[ c(28:32) ]
meta_data_2$sample<-rownames(meta_data_2)
diversity <- merge.data.frame(observed.diversity, shannon.diversity, by='sample')
rownames(diversity)<-diversity$sample
meta_diversity <- merge.data.frame(meta_data_2, diversity, by='sample')

kruskal.test(observed_features ~ case_2, data = meta_diversity)
kruskal.test(shannon_entropy ~ case_2, data = meta_diversity)

div_pos <- meta_diversity[meta_diversity$case_2=='Positive',]
div_neg <- meta_diversity[meta_diversity$case_2=='Negative',]

quantile(div_pos$observed_features)
quantile(div_neg$observed_features)

quantile(div_pos$shannon_entropy)
quantile(div_neg$shannon_entropy)

observed.diversity$measure <-"Observed features"
names(observed.diversity)<-c('value', 'sample', 'measure')
shannon.diversity$measure <- "Shannon entropy"
names(shannon.diversity)<-c('value', 'sample', 'measure')

diversity <- rbind(observed.diversity, shannon.diversity)
diversity$case <- meta_data_2$case_2[match(diversity$sample, meta_data_2$sample)]
diversity$case <- factor(diversity$case, levels = c('Negative', 'Positive'), labels = c('uninfected', 'infected'))
diversity$measure <- factor(diversity$measure, levels = c('Shannon entropy', 'Observed features'))

s

ggplot(data = diversity, aes(x=case, y=value, fill=case))+
  geom_boxplot() + facet_wrap(~measure, scale='free')+
  theme(legend.position = 'none')+
  xlab("S. mansoni infection status")

##PANEL B and C
asvs_abundance <- read.table("~/Novogene/pet-table/feature-tablew-clean.tsv", header = TRUE, row.names=1, comment.char = "#")
with_diet_laviiswa<-read.csv("~/Novogene/pet-table/allcombineddata.csv")
#asvs_abundance <- read.table("feature-table.tsv", header = TRUE, row.names=1)
#asvs_abundance <- import_biom(BIOMfilename = "feature-table.biom", treefilename = "tree.nwk")
#asvs_abundance <- read.table("feature-table.tsv", header = TRUE, row.names=1)
asv_taxonomy <- read.delim("~/Novogene/pet-table/taxonomy.tsv", row.names = 1)
asv_taxonomy[] <- lapply(asv_taxonomy, function(x) gsub("d__|p__|c__|o__|f__|g__|s__|\\[|\\]","", x))
tmp <- cSplit(asv_taxonomy, "Taxon", ";")
tmp$Confidence<-NULL
names(tmp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")  #assign new column names
colnames(meta_data)[colnames(meta_data) == "case_2"] ="Hel_inf_status"
tax <- data.frame(tmp, row.names = rownames(asv_taxonomy))
#meta_data <- read.delim("~/Novogene/With_case_final_novogene_received_with_lab_ids_and_outcomedata.txt", header=TRUE)
rownames(meta_data)<-meta_data$sample.id
snames<-colnames(asvs_abundance)[colnames(asvs_abundance)%in%rownames(meta_data)]
asvs_abundance<-asvs_abundance[snames]
meta_data<-meta_data[rownames(meta_data)%in%snames,]
ASV <- otu_table(as.matrix(asvs_abundance), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
meta_data$sample <- row.names(meta_data)
colnames(meta_data)[colnames(meta_data) == "case_2"] ="Hel_inf_status"
colnames(meta_data)[colnames(meta_data) == "cholesterol2"] ="Chol_T"
colnames(meta_data)[colnames(meta_data) == "ldl_c"] ="Chol_L"
colnames(meta_data)[colnames(meta_data) == "dbpnew"] ="BP_Dia"
colnames(meta_data)[colnames(meta_data) == "sbpnew"] ="BP_Sys"
SAM <- sample_data(meta_data)
physeq <- merge_phyloseq(phyloseq(ASV, TAX), SAM)
physeq
#saveRDS(object = physeq, file = 'PhyloseqObjectHelminths.rds')
######
##pickouttheruralguys
Ruralphyseq<-subset_samples(physeq, !is.na(sample_data(physeq)$Rural_infected))
Urbanphyseq<-subset_samples(physeq, !is.na(sample_data(physeq)$Urban_infected))

###specificgroups
ruralphyseqcase<-subset_samples(Ruralphyseq, !is.na(sample_data(Ruralphyseq)$Hel_inf_status))
PCoA.ord.bray <- ordinate(ruralphyseqcase, "PCoA", "bray")
#NMDS.ord.bray <- ordinate(ruralphyseqcase, "NMDS", "bray")
sample_data(ruralphyseqcase)$Hel_inf_status<-factor(sample_data(ruralphyseqcase)$Hel_inf_status,
                                                    levels = c('Negative', 'Positive'),
                                                    labels = c('uninfected','infected'))
PCoA.ord <- plot_ordination(ruralphyseqcase, PCoA.ord.bray, color = "Hel_inf_status")
#PCoA.ord <- plot_ordination(ruralphyseqcase, NMDS.ord.bray, color = "Hel_inf_status")
PCoA.ord <- PCoA.ord + theme(axis.text = element_text(size = 16, face = "bold"),axis.title = element_text(size = 18, face = "bold"), legend.title = element_text(size = 14)) +
  theme_bw() + labs(color = "Hel_inf_status") + geom_point(size = 5)+stat_ellipse()+ ggtitle("Sm+ve vs Sm-ve Beta diversity dissimilarity")
PCoA.ord$labels$colour<-'S. mansoni infectious status'
PCoA.ord

bray_dist<-phyloseq::distance((ruralphyseqcase), method = "bray")
tmp<-data.frame(sample_data(ruralphyseqcase))
adonis(bray_dist~Hel_inf_status, data=tmp)

###quantify separation
require(vegan)

###specificgroups
Urbanphyseqcase<-subset_samples(Urbanphyseq, !is.na(sample_data(Urbanphyseq)$Hel_inf_status))
PCoA.ord.bray <- ordinate(Urbanphyseqcase, "PCoA", "bray")
sample_data(Urbanphyseqcase)$Hel_inf_status<-factor(sample_data(Urbanphyseqcase)$Hel_inf_status,
                                                    levels = c('Negative', 'Positive'),
                                                    labels = c('uninfected','infected'))
PCoA.ord <- plot_ordination(Urbanphyseqcase, PCoA.ord.bray, color = "Hel_inf_status")
PCoA.ord <- PCoA.ord + theme(axis.text = element_text(size = 16, face = "bold"),axis.title = element_text(size = 18, face = "bold"), legend.title = element_text(size = 14)) +
  theme_bw() + labs(color = "Hel_inf_status") + geom_point(size = 5)+stat_ellipse()+ ggtitle("Sm+ve vs Sm-ve Beta diversity dissimilarity")
PCoA.ord$labels$colour<-'S. mansoni infectious status'
PCoA.ord

bray_dist<-phyloseq::distance((Urbanphyseq), method = "bray")
tmp<-data.frame(sample_data(Urbanphyseq))
adonis(bray_dist~Hel_inf_status, data=tmp)


#### PANEL D AND E

# Rural 
counts_rural <- as.matrix(otu_table(Ruralphyseq))
colDataRural <- sample_data(Ruralphyseq)
colnames(colDataRural)[colnames(colDataRural) == "case_2"] ="Hel_inf_status"

se_rural<-SummarizedExperiment(assays=list(counts_rural),colData=colDataRural)
se_rural

res_rural_infected <- lefser(se_rural, groupCol = "Rural_infected", lda.threshold=2.0)
res_rural_infected$Names<-tax$Genus[match(res_rural_infected$Names, rownames(tax))]
lefserPlot(res_rural_infected)

# Urban 
counts_urban <- as.matrix(otu_table(Urbanphyseq))
colDataUrban <- sample_data(Urbanphyseq)
colnames(colDataRural)[colnames(colDataRural) == "case_2"] ="Hel_inf_status"

se_urban<-SummarizedExperiment(assays=list(counts_urban),colData=colDataUrban)
se_urban

res_urban_infected <- lefser(se_urban, groupCol = "Urban_infected", lda.threshold=2.0)
res_urban_infected$Names<-tax$Genus[match(res_urban_infected$Names, rownames(tax))]
lefserPlot(res_urban_infected)
####
merged_data <- merge(
  subset_data,
  meta_data,
  by = "slabno",
  all = FALSE   # keeps only matching slabno in both
)

# Check result
dim(merged_data)
head(merged_data)
