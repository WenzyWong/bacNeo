##########################################
# Under the project: bacNeo
# Running test for AEG-peptides
#
# Yunzhe Wang, yunzhewang24@m.fudan.edu.cn
# Created: 2025-06-30
# Updated: 2025-07-05
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
library(forcats)

# Data preprocessing
pep.raw <- readRDS("./rds/peptides_raw.rds")
cpm.species <- readRDS("./rds/rna_cpm.rds")
strong.binders <- readRDS("./rds/strong_binder.rds")
hla.samples <- readRDS("./rds/patient_hla.rds")

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

# Plotting HLA allele distribution
draw_allele <- hla.samples[ , c(2:ncol(hla.samples))] %>%
  rename(allele = Type) %>%
  mutate(allele = gsub("N", "", allele),
         allele_frequency = Count / 188) %>%
  mutate(allele_frequency = 100 * allele_frequency) %>%
  arrange(Group, allele_frequency)

draw_allele <- draw_allele[!duplicated(draw_allele$allele),]
p_allele_distribution <-
  ggplot(draw_allele, aes(x = allele_frequency,
                        y = factor(allele, levels = unique(allele)), 
                        fill = Group)) + 
  facet_wrap(~Group, scales = "free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#B574AE", "#F8C77A", "#6387C5")) +
  ylab("HLA alleles") + 
  xlab("Percentage (%)") + 
  geom_vline(xintercept = 10, colour = "red3") + 
  annotate("text", x = 10, y = 2, label = "10%", colour = "red3") +
  theme_bw() +
  theme(axis.text = element_text(colour = 1))
pdf("./outputs/allele_distribution.pdf", width = 10, height = 7)
print(p_allele_distribution)
dev.off()

# Plotting scatters
allele_freq <- draw_allele %>%
  group_by(allele) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / nrow(draw_allele) * 100)

total_peptides <- n_distinct(strong.binders$Peptide)

peptide_binding <- strong.binders %>%
  group_by(Allele) %>%
  summarise(peptide_count = n_distinct(Peptide)) %>%
  mutate(peptide_percentage = peptide_count / total_peptides * 100)

allele_peptide_data <- allele_freq %>%
  inner_join(peptide_binding, by = c("allele" = "Allele")) %>%
  arrange(desc(peptide_percentage))

rainbow_colour <- paletteer_d("khroma::smoothrainbow")
n_allele <- length(unique(allele_peptide_data$allele))
if (n_allele < length(rainbow_colour)) {
  allele_colour <- rainbow_colour[seq(1, length(rainbow_colour), by = round(length(rainbow_colour)/n_allele))]
} else {
  all_colour <- paletteer_d("palettesForR::Named")
  allele_colour <- all_colour[seq(1, length(all_colour), by = round(length(all_colour)/n_allele)-1)]
}

p_scatter <-
  ggplot(allele_peptide_data, aes(x = percentage, y = peptide_percentage,
                                  color = allele)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_colour_manual(values = allele_colour) + 
  ggrepel::geom_label_repel(aes(label = allele),
                            nudge_x = 1,
                            nudge_y = 1,
                            alpha = .8) +
  labs(
    x = "HLA alleles (%)",
    y = "Bound peptides (%)"
  ) +
  theme_test() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
pdf("./outputs/scatter_allele_peptide_percentage.pdf", width = 8, height = 5)
print(p_scatter)
dev.off()

# Plotting Jaccard index
p_jaccard <-
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
print(p_jaccard)
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

p_heat <-
  Heatmap(sign(pep_detected[high_detection_species, ]),
          name = "Detection", col = col_fun, 
          column_title = "Samples", column_title_side = "bottom",
          show_column_names = F,
          row_names_gp = gpar(fontface = "italic", fontsize = 10)
  )
pdf("./outputs/detection_matrix.pdf", width = 6, height = 4)
draw(p_heat)
dev.off()

# Peptide selection
strong.binders <- strong.binders %>%
  arrange(BArank) %>%
  add_count(Allele, name = "AlleleFreq")
write.csv(strong.binders, "./outputs/strong_binder_stats.csv")
table(strong.binders$Allele)
binder.sankey <- data.frame(
  Source = c(strong.binders$Allele, strong.binders$Peptide, strong.binders$ProteinID),
  Target = c(strong.binders$Peptide, strong.binders$ProteinID, strong.binders$Species)
)

binder.sankey <- binder.sankey %>%
  mutate(Value = c(table(Source)[Source]))
dim(binder.sankey)
write.csv(binder.sankey, "./outputs/strong_binders_sankey.csv")

# Plot chord diagram
binder.adjacency <- binder.sankey
colnames(binder.adjacency) <- c("from", "to", "value")

node_values <- aggregate(value ~ from, data = binder.adjacency, sum)
names(node_values) <- c("node", "total_value")
node_values2 <- aggregate(value ~ to, data = binder.adjacency, sum)
names(node_values2) <- c("node", "total_value")
all_nodes <- rbind(node_values, node_values2)
node_totals <- aggregate(total_value ~ node, data = all_nodes, sum)

show_labels <- node_totals$node[node_totals$total_value > 10]

set.seed(42)
colors <- sample(rainbow(n_sectors))
names(colors) <- all_sectors

pdf("./outputs/chord_diagram.pdf", width = 10, height = 10)
circos.par(start.degree = 180)
chordDiagram(binder.adjacency, 
             grid.col = colors,
             annotationTrack = "grid",
             annotationTrackHeight = 0.05)

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% show_labels) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.2, 
                sector.name, facing = "clockwise", 
                niceFacing = TRUE, adj = c(0, 0.5),
                cex = 0.8)
  }
}, bg.border = NA)

circos.clear()
dev.off()
