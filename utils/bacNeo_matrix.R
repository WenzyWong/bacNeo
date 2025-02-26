#!/usr/bin/env Rscript
#################################################################################
# Rscript for combining bacterial counts / CPM / abundance of different samples
# and creating a matrix

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-02-26
#################################################################################
if (!requireNamespace("argparse", quietly = TRUE)) {
    install.packages("argparse", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    remotes::install_github("corybrunson/ggalluvial@main", build_vignettes = TRUE)
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork", repos = "https://cloud.r-project.org")
}

library(argparse)
library(ggplot2)
library(ggalluvial)
library(patchwork)

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
    help = "Data normalization method (e.g., 'CPM', 'abundance')",
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
            file = file.path(DIR_RES, paste0("matrix_", NORM, "_", TAXONOMY, ".txt"))
            )

##############################
# Visualization
if (NORM == "abundance" | NORM == "CPM") {
  mtx <- read.delim(file.path(DIR_RES, paste0("matrix_", NORM, "_", TAXONOMY, ".txt")), row.names = 1)
  if (NORM == "CPM") {
    abund_select <- sort(apply(mtx / (1e+4), MARGIN = 1, FUN = mean), decreasing = T)[1:10]
  } else {
    abund_select <- sort(apply(mtx, MARGIN = 1, FUN = mean), decreasing = T)[1:10]
  }
  abund_select <- data.frame(
    taxonomy = c(names(abund_select), "others"),
    abundance = c(abund_select, "others" = 100 - sum(abund_select)))

  abund_alluvium <- mtx[abund_select$taxonomy[1:10], ]
  abund_alluvium <- rbind(
    abund_alluvium,
    "others" = 100 - apply(abund_alluvium, MARGIN = 2, FUN = sum)
  )
  abund_alluvium$taxonomy <- rownames(abund_alluvium)
  abund_alluvium <- reshape2::melt(abund_alluvium)

  abund_col <- c("#64A4CC", "#9CCEE3", "#ADE0DC", "#D3F2D6", "#ECF6C8", 
                "#FEEDAA", "#F89D59", "#E75B3A", "#CD2626", "#9A342C", 
                "#DCDCDC")

  p1 <- 
    ggplot(data = abund_alluvium, aes(x = variable, y = value, 
                                    alluvium = factor(taxonomy, 
                                                      levels = unique(taxonomy)))) +
    geom_alluvium(aes(fill = factor(taxonomy, 
                                    levels = unique(taxonomy))), 
                  alpha = 1) +
    scale_fill_manual(values = abund_col, name = "Taxonomy") +
    xlab("Samples") +
    ylab("Relative abundance (%)") +
    theme_test() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, colour = 1), 
          axis.text.y = element_text(colour = 1),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="none")

  p2 <-
    ggplot(data = abund_select, aes(x = 1, y = abundance,
                                  fill = factor(taxonomy, levels = taxonomy))) + 
    geom_bar(width = .5, stat = "identity", position = "stack") +
    scale_fill_manual(values = abund_col, name = "Taxonomy") +
    xlab("Sample (average)") +
    ylab("Relative abundance (%)") +
    theme_test() +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(colour = 1),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

  p <- p1 + p2 + plot_layout(guides = "collect", widths = c(6, 2))
  pdf(file.path(DIR_RES, paste0("Plot_abundance_", TAXONOMY, ".pdf")), width = 6, height = 6)
  print(p)
  dev.off()
}
