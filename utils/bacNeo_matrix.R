#!/usr/bin/env Rscript
#################################################################################
# Rscript for combining bacterial counts / CPM / abundance of different samples
# and creating a matrix

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-01-12
#################################################################################
if (!requireNamespace("argparse", quietly = TRUE)) {
    install.packages("argparse", repos = "https://cloud.r-project.org")
}

library(argparse)

###################
# Parsing arguments
parser <- ArgumentParser(
    description = "Process and merge species data from multiple samples.",
    formatter_class = "argparse.RawTextHelpFormatter"
)
# Add arguments
parser$add_argument("-d", "--dir", 
    type = "character",
    help = "Root directory path for result files",
    metavar = "DIR_PATH"
)

parser$add_argument("-l", "--level", 
    type = "character",
    help = "Taxonomic level (e.g., 's' for species level)",
    metavar = "TAXONOMY"
)

parser$add_argument("-n", "--norm", 
    type = "character",
    help = "Data normalization method (e.g., 'CPM')",
    metavar = "NORM"
)

# Parse command line arguments
args <- parser$parse_args()
if (is.null(args$dir) || is.null(args$level) || is.null(args$norm)) {
    parser$print_help()
    stop("All required parameters must be provided. Use -h to see help.", call. = FALSE)
}

DIR_RES <- args$dir
TAXONOMY <- args$level
NORM <- args$norm

######
# Main
SAMPLE <- list.files(DIR_RES, recursive = F)

dir_norm <- paste0(DIR_RES, "/", SAMPLE, "/normalized_",
                   TAXONOMY, ".txt")

result <- NULL

for (i in seq_along(dir_norm)) {
  if (!file.exists(dir_norm[i])) {
    warning(sprintf("File does not exist %s", dir_norm[i]))
    next
  }
  
  curr_data <- read.delim(dir_norm[i])
  
  curr_data <- curr_data[, c("species", NORM)]
  colnames(curr_data)[2] <- SAMPLE[i]
  
  if (is.null(result)) {
    result <- curr_data
  } else {
    result <- merge(result, curr_data, by = "species", all = T)
  }
}
result[is.na(result)] <- 0

result <- result[order(rowMeans(result[, -1]), decreasing = T), ]

write.table(result, sep = "\t", row.names = F, quote = F,
            file = paste0(DIR_RES, "/matrix_", NORM, "_", TAXONOMY, ".txt")
            )