##########################################
# Under the project: bacNeo
# Running test for scRNA-seq in GC
#
# Yunzhe Wang, yunzhewang24@m.fudan.edu.cn
# Created: 2025-07-21
# Updated: 2025-07-31
##########################################
library(Seurat)
library(ggplot2)
library(dplyr)

meta <- read.csv("./SraRunTable.csv")
rownames(meta) <- meta$Run
id.sample <- SraRunTable$Library.Name

filt.norm <- readRDS("./rds/scRNA_filtered_normed.rds")

col.celltype <- c("Parietal" = "#6a3d9a",
                  "Surface mucous" = "#9B62A7FF",
                  "Fibroblast" = "#C3A8D1FF",
                  "B" = "#4E79C5FF",
                  "Plasma" = "#4E96BCFF",
                  "APC" = "#8dd3c7",
                  "T" = "#01665e",
                  "Endothelial" = "#b2df8a",
                  "Epithelial" = "#fdbf6f",
                  "Enteroendocrine" = "#ff7e00",
                  "Stem-like" = "tomato",
                  "Malignant epithelial" = "#DF4828FF",
                  "Proliferate" = "#95211BFF"
)
DimPlot(filt.norm, reduction = "tsne", group.by = c("cell.type"), cols = col.celltype, label = T)
# 01_overall_celltypes.pdf, 6.5 * 5

# Load kraken outputs
# Change the path down below to your bacc output directory
kk.files <- list.files("./", pattern = "KRAKEN", recursive = T)
for (i in 1:length(kk.files)) {
  tmp.file <- kk.files[i]
  tmp.sam <- meta[gsub("/..*", "", tmp.file), "Library.Name"]
  tmp.out <- read.delim(paste0("./species/", tmp.file), header = F) %>%
    filter(V3 != "9606") %>%
    mutate(barcode = sapply(strsplit(V2, "\\."), function(x) x[3])) %>%
    mutate(row = paste0(tmp.sam, "_", barcode, "-1"))
  if (i == 1) {kk.out <- tmp.out} else {kk.out <- rbind(kk.out, tmp.out)}
}

# Extract bacteria (ncbi taxid == 2)
bac.all <- kk.out %>%
  filter(V3 == "2")
dim(bac.all)
bac.all.ct <- bac.all %>%
  group_by(row) %>%
  summarise(total_count = n(), .groups = "drop")
dim(bac.all.ct)

match_indices <- match(rownames(filt.norm@meta.data), bac.all.ct$row)
length(match_indices)

filt.norm@meta.data$bac_all <- ifelse(
  is.na(match_indices), 0, log10(bac.all.ct$total_count[match_indices])
)
FeaturePlot(filt.norm, reduction = "tsne", features = c("bac_all")) +
  scale_colour_distiller(palette = "Blues", direction = 0)
