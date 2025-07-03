##########################################
# Under the project: bacNeo
# Running test for AEG-peptides
#
# Yunzhe Wang, yunzhewang24@m.fudan.edu.cn
# Created: 2025-06-30
# Updated: 2025-07-02
##########################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(corrplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(VennDiagram)
library(reshape2)
library(scales)

# Data preprocessing
pep.raw <- readRDS("./rds/peptides_raw.rds")
cpm.species <- readRDS("./rds/rna_cpm.rds")
strong.binders <- readRDS("./rds/strong_binder.rds")

mtx.pep <- sign(pep.raw[ , 4:(ncol(pep.raw) - 1)])
rownames(mtx.pep) <- pep.raw$Protein.IDs

pep.lands <- mtx.pep %>%
  mutate(species = pep.raw[rownames(mtx.pep), "Species"])
pep.agg <- aggregate(pep.lands[1:ncol(mtx.pep)],
                     by = list(pep.lands$species),
                     FUN = sum)
rownames(pep.agg) <- pep.agg$Group.1
pep.agg <- as.matrix(pep.agg[ , -1])

# Extracting overlaps
common_samples <- intersect(colnames(cpm.species), colnames(pep.agg))
common_species <- intersect(rownames(cpm.species), rownames(pep.agg))
cpm_subset <- cpm.species[common_species, common_samples]
pep_subset <- pep.agg[common_species, common_samples]

# Filtering
rna_threshold <- 10
peptide_threshold <- 1

rna_detected <- cpm_subset > rna_threshold
pep_detected <- pep_subset > peptide_threshold

# Stats
detection_stats <- data.frame(
  species = common_species,
  rna_detection_freq = rowSums(rna_detected, na.rm = TRUE) / ncol(rna_detected),
  pep_detection_freq = rowSums(pep_detected, na.rm = TRUE) / ncol(pep_detected),
  both_detected_freq = rowSums(rna_detected & pep_detected, na.rm = TRUE) / ncol(rna_detected),
  either_detected_freq = rowSums(rna_detected | pep_detected, na.rm = TRUE) / ncol(rna_detected),
  stringsAsFactors = FALSE
)

detection_stats$jaccard_index <- detection_stats$both_detected_freq / detection_stats$either_detected_freq
detection_stats$detection_concordance <- detection_stats$both_detected_freq / 
  pmax(detection_stats$rna_detection_freq, detection_stats$pep_detection_freq)

# Plotting Jaccard index
p1 <-
  ggplot(detection_stats, aes(x = rna_detection_freq, y = pep_detection_freq)) +
  geom_point(aes(size = both_detected_freq, color = jaccard_index), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_viridis_c(name = "Jaccard\nIndex") +
  scale_size_continuous(range = c(1, 8), name = "Co-detection\nFrequency") +
  labs(
    title = "Detection Consistency",
    x = "RNA Detection Frequency",
    y = "Peptide Detection Frequency",
  ) +
  theme_bw() +
  geom_text_repel(
    data = detection_stats[detection_stats$jaccard_index > 0.3, ],
    aes(label = species), size = 3, max.overlaps = 15
  ) +
  theme(axis.text = element_text(colour = 1))
pdf("./outputs/detection_consistency.pdf", width = 6.2, height = 5)
print(p1)
dev.off()

# Plotting species detected consistently
high_detection_species <- detection_stats$species[
  detection_stats$rna_detection_freq > 0.2 & detection_stats$pep_detection_freq > 0.2
  ]

detection_matrix <- rbind(
  RNA = ifelse(rna_detected[high_detection_species, ], 1, 0),
  Peptide = ifelse(pep_detected[high_detection_species, ], 1, 0)
)

sample_variance <- apply(detection_matrix, 2, function(x) var(x, na.rm = TRUE))
valid_samples <- names(sample_variance)[!is.na(sample_variance) & sample_variance > 0]

print(length(high_detection_species))
print(length(valid_samples))

detection_matrix_filtered <- detection_matrix[, valid_samples]

col_fun <- c("0" = "white", "1" = "darkblue")

p2 <-
  Heatmap(sign(pep_detected[high_detection_species, ]),
          name = "Detection", col = col_fun, 
          column_title = "Samples", column_title_side = "bottom",
          show_column_names = F,
          row_names_gp = gpar(fontface = "italic", fontsize = 10)
  )
pdf("./outputs/detection_matrix.pdf", width = 6, height = 4)
draw(p2)
dev.off()

# Peptide selection
table(strong.binders$Allele)
binder.sankey <- data.frame(
  Source = c(strong.binders$Allele, strong.binders$Peptide, strong.binders$ProteinID),
  Target = c(strong.binders$Peptide, strong.binders$ProteinID, strong.binders$Species)
)

binder.sankey <- binder.sankey %>%
  mutate(Value = c(table(Source)[Source]))
dim(binder.sankey)
write.csv(binder.sankey, "./outputs/strong_binders_sankey.csv")
