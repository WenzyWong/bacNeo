#################################################################################
# Rscript for calculating CPM and abundance of bacterial reads count

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-01-11
#################################################################################
args = commandArgs(trailingOnly=TRUE)

DIR_RES <- args[1]
TAXONOMY <- args[2]

library(edgeR)

dir_counts <- paste0(DIR_RES, "/", "counts_", TAXONOMY, ".txt")
dt_counts <- read.delim(dir_counts, header = F)
colnames(dt_counts) <- c("taxa", "counts")

counts_matrix <- matrix(dt_counts$counts, ncol = 1)
rownames(counts_matrix) <- dt_counts$taxa

dge <- DGEList(counts = counts_matrix)

# Calculating CPM
cpm_norm <- cpm(dge)

# Calculating abundance
abundance <- cpm_norm / 1e4

results <- data.frame(
  taxa = dt_counts$taxa,
  raw_counts = dt_counts$counts,
  abundance = as.numeric(abundance),
  CPM = as.numeric(cpm_norm)
)

write.table(results,
            file = paste0(DIR_RES, "/long_norm_", TAXONOMY, ".txt"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
