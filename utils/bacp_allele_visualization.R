setwd("/data/yzwang/tests/visualization_250226/")
DIR <- "../neo_merge_250119/Res/bacp/00_allele_summary/"
FILE <- list.files(DIR)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("forcats", quietly = TRUE)) {
  install.packages("forcats", repos = "https://cloud.r-project.org")
}

library(dplyr)
library(ggplot2)

allele_dt <- data.frame()
for (i in 1:length(FILE)){
  tmp_path <- file.path(DIR, FILE[i])
  tmp_dt <- read.table(tmp_path)
  allele_dt <- rbind(allele_dt, tmp_dt)
  rownames(allele_dt)[(2 * i - 1):(2 * i)] <- paste0(gsub(".txt", "", FILE[i]), "_", 1:2)
}

draw_allele <- allele_dt %>%
  rename(allele = V1) %>%
  mutate(sample = gsub("_..*", "", rownames(.)))

draw_line <- 0.2 * nrow(draw_allele)

ggplot(draw_allele, aes(x = forcats::fct_infreq(allele))) + 
  geom_bar(fill = "grey") +
  xlab("HLA alleles") + 
  ylab("Count") + 
  geom_hline(yintercept = draw_line, colour = "red3") + 
  annotate("text", x = .6, y = draw_line, label = "20%", colour = "red3")
  theme_bw() +
  theme(axis.text = element_text(colour = 1),
        axis.text.x = element_text(angle = 90, vjust = .5))

