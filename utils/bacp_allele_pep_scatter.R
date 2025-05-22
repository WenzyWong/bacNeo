args = commandArgs(trailingOnly=TRUE)

DIR <- args[1]

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("paletteer", quietly = TRUE)) {
  install.packages("paletteer", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos = "https://cloud.r-project.org")
}

library(dplyr)
library(ggplot2)
library(paletteer)

strong <- read.csv(file.path(DIR, "Strong_binders.csv"))
draw_allele <- read.csv(file.path(DIR, "00_allele_summary", "Allele_summary.csv"))

allele_freq <- draw_allele %>%
  group_by(allele) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / nrow(draw_allele) * 100)

total_peptides <- n_distinct(strong$Peptide_ID)

peptide_binding <- strong %>%
  group_by(HLA_allele) %>%
  summarise(peptide_count = n_distinct(Peptide_ID)) %>%
  mutate(peptide_percentage = peptide_count / total_peptides * 100)

allele_peptide_data <- allele_freq %>%
  inner_join(peptide_binding, by = c("allele" = "HLA_allele")) %>%
  arrange(desc(peptide_percentage))

rainbow_colour <- paletteer_d("khroma::smoothrainbow")
n_allele <- length(unique(allele_peptide_data$allele))
if (n_allele < length(rainbow_colour)) {
  allele_colour <- rainbow_colour[seq(1, length(rainbow_colour), by = round(length(rainbow_colour)/n_allele))]
} else {
  all_colour <- paletteer_d("palettesForR::Named")
  allele_colour <- all_colour[seq(1, length(all_colour), by = round(length(all_colour)/n_allele))]
}

p <- ggplot(allele_peptide_data, aes(x = percentage, y = peptide_percentage,
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
pdf(file.path(DIR, "Scatter_allele_peptide_percentage.pdf"), width = 6, height = 5)
print(p)
dev.off()