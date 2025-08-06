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
evenness.diversity <- read.delim("~/Novogene/core-metrics-results/evenness_vector.tsv/alpha-diversity.tsv", row.names = 1)
evenness.diversity$sample<-rownames(evenness.diversity)

meta_data <- read.delim("~/Novogene/With_case_final_novogene_received_with_lab_ids_and_outcomedata.txt", header=TRUE)
rownames(meta_data)<-meta_data$sample.id
meta_data_2 <- meta_data
meta_data_2$sample<-rownames(meta_data_2)
meta_data_2 <- meta_data_2[c('setting.x','sample')]
diversity <- merge.data.frame(observed.diversity, shannon.diversity, by='sample')
diversity <- merge.data.frame(diversity, evenness.diversity, by='sample')

rownames(diversity)<-diversity$sample
meta_diversity <- merge.data.frame(meta_data_2, diversity, by='sample')

kruskal.test(observed_features ~ setting.x, data = meta_diversity)
kruskal.test(shannon_entropy ~ setting.x, data = meta_diversity)
kruskal.test(pielou_evenness ~ setting.x, data = meta_diversity)

observed.diversity$measure <-"Observed features"
names(observed.diversity)<-c('value', 'sample', 'measure')
shannon.diversity$measure <- "Shannon entropy"
names(shannon.diversity)<-c('value', 'sample', 'measure')
evenness.diversity$measure <- "Pielou evenness"
names(evenness.diversity)<-c('value', 'sample', 'measure')

diversity <- rbind(observed.diversity, shannon.diversity,evenness.diversity)
diversity$setting <- meta_data_2$setting.x[match(diversity$sample, meta_data_2$sample)]
diversity$setting <- factor(diversity$setting, levels = c('Rural', 'Urban'), labels = c('Rural', 'Urban'))
diversity$measure <- factor(diversity$measure, levels = c('Shannon entropy', 'Observed features', 'Pielou evenness'))

ggplot(data = diversity, aes(x=setting, y=value, fill=setting))+
  geom_boxplot() + facet_wrap(~measure, scale='free')+
  theme(legend.position = 'none')+
  xlab("Setting")

###### PANEL B ########

ruralurabnphyseq <- readRDS('/Users/mbmhvbw2/Novogene/codeNaturecomms/ruralurabnphyseq.rds')
PCoA.ord.bray <- ordinate(ruralurabnphyseq, "PCoA", "bray")
#NMDS.ord.bray <- ordinate(ruralurabnphyseq, "NMDS", "bray")
PCoA.ord <- plot_ordination(ruralurabnphyseq, PCoA.ord.bray, color = "setting.x")
#PCoA.ord <- plot_ordination(ruralurabnphyseq, NMDS.ord.bray, color = "setting.x")
PCoA.ord <- PCoA.ord + theme(axis.text = element_text(size = 16, face = "bold"),axis.title = element_text(size = 18, face = "bold"), legend.title = element_text(size = 14)) +
  theme_bw() + labs(color = "setting.x") + geom_point(size = 2)+stat_ellipse()+ ggtitle("Rural Vs Urban Beta-diversity dissimilarity") 
PCoA.ord
###quantify separation
bray_dist<-phyloseq::distance((ruralurabnphyseq), method = "bray")
tmp<-data.frame(sample_data(ruralurabnphyseq))
adonis2(bray_dist~setting.x, data=tmp)

###Panel C #####

setwd("/Users/mbmhvbw2/Novogene/pet-table")
#######
# Load required libraries
library(SummarizedExperiment)
library(GenomicRanges)
library(lefser)
library(splitstackshape)
library(ggplot2)

# Step 1: Load and process metadata
meta_data <- read.delim("~/Novogene/With_case_final_novogene_received_with_lab_ids_and_outcomedata.txt", header = TRUE)
rownames(meta_data) <- meta_data$sample.id
#colnames(meta_data)[colnames(meta_data) == "case_2"] <- "Hel_inf_status"  # if needed
colData <- DataFrame(meta_data)

# Step 2: Load ASV abundance table and align with metadata
asvs_abundance <- read.table("feature-tablew-clean.tsv", header = TRUE, row.names = 1, comment.char = "#")
asvs_abundance <- asvs_abundance[, rownames(meta_data)]  # ensure sample alignment

# Step 3: Convert to SummarizedExperiment
counts <- as.matrix(unname(asvs_abundance))
se <- SummarizedExperiment(assays = list(counts = counts), colData = colData)

# Step 4: Load and clean taxonomy
asv_taxonomy <- read.delim("taxonomy.tsv", row.names = 1, stringsAsFactors = FALSE)
asv_taxonomy[] <- lapply(asv_taxonomy, function(x) gsub("d__|p__|c__|o__|f__|g__|s__|\\[|\\]", "", x))
tmp <- cSplit(asv_taxonomy, "Taxon", ";", type.convert = FALSE)
tmp$Confidence <- NULL
colnames(tmp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- data.frame(tmp, row.names = rownames(asv_taxonomy))

# Step 5: Run LEfSe
# Change 'setting.x' to the actual grouping variable name in your metadata
if (!"setting.x" %in% colnames(meta_data)) stop("Column 'setting.x' not found in metadata.")
res <- lefser(se, groupCol = "setting.x", lda.threshold = 2.0)

# Step 6: Replace ASV IDs with genus-level names (or fallback to ASV ID if Genus is NA)
res$Genus <- tax$Genus[match(res$Names, rownames(tax))]
res$Genus[is.na(res$Genus) | res$Genus == ""] <- res$Names[is.na(res$Genus) | res$Genus == ""]
resF<-res
res$Genus <- factor(res$Genus, levels = res$Genus[order(res$LDA, decreasing = TRUE)])
# Reorder Genus levels by scores (which are your LDA values)
res$Genus <- factor(res$Genus, levels = res$Genus[order(res$scores, decreasing = TRUE)])
res$Genus_unique <- make.unique(res$Genus)
res$Genus_unique <- factor(res$Genus_unique, levels = res$Genus_unique[order(res$scores, decreasing = TRUE)])
####
library(dplyr)
library(ggplot2)
library(stringr)

# Remove unnamed taxa (i.e., Genus names that look like long hashes)
res_clean <- res %>%
  filter(!grepl("^[0-9a-f]{32}$", Genus))  # excludes 32-character hex strings

# Keep only one representative per Genus (highest |LDA score|)
res_unique <- res_clean %>%
  group_by(Genus) %>%
  slice_max(order_by = abs(scores), n = 1, with_ties = FALSE) %>%
  ungroup()

# Wrap long names
res_unique$Genus_wrapped <- stringr::str_wrap(res_unique$Genus, width = 25)

# Set factor order
res_unique$Genus_wrapped <- factor(res_unique$Genus_wrapped, levels = res_unique$Genus_wrapped[order(res_unique$scores)])

# Plot
ggplot(res_unique, aes(x = Genus_wrapped, y = scores, fill = scores > 0)) +
  geom_bar(stat = "identity", width = 0.75, color = "black") +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#33a02c", "FALSE" = "#e31a1c"),  # Blue = Enriched, Green = Depleted
    labels = c("Enriched", "Depleted")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold.italic", size = 9, margin = margin(r = 5)),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  labs(
    title = "LEfSe: Differentially Abundant Genera Between Groups",
    x = "Genus",
    y = "LDA Score (log10)"
  )

###

