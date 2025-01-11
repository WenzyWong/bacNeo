#################################################################################
# Rscript for calculating CPM and abundance of bacterial reads count

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-01-11
#################################################################################
args = commandArgs(trailingOnly=TRUE)

DIR_RES <- args[1]
SAMPLE <- args[2]
TAXONOMY <- args[3]

library(edgeR)

dir_counts <- paste0(DIR_RES, "/", samples[i], "/", "counts_", TAXONOMY, ".txt")
dt_counts <- read.delim(tmp_counts_path, header = F)
colnames(dt_counts) <- c("species", "counts")

counts_matrix <- matrix(dt_counts$counts, ncol = 1)
rownames(counts_matrix) <- dt_counts$species
colnames(counts_matrix) <- SAMPLE

dge <- DGEList(counts = counts_matrix)

# Calculating CPM
cpm_norm <- cpm(dge)

# Calculating abundance
abundance <- cpm_values / 1e4

results <- data.frame(
  species = dt_counts$species,
  raw_counts = dt_counts$counts,
  abundance = as.numeric(abundance),
  CPM = as.numeric(cpm_norm)
)

write.table(results, 
            file = paste0(DIR_RES, "/", SAMPLE, "/normalized_", TAXONOMY, ".txt"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
