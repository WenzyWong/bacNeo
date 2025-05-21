#################################################################################
# Rscript for extracting strong and weak binders from netMHCpan output files
# and calculating logIC50 for predicted TAP-binding efficiency

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-05-16
#################################################################################
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
OUTPUT <- args[1]

hla_files <- list.files(path = OUTPUT, pattern = "HLA-", recursive = T)

all_strong <- data.frame()
all_weak <- data.frame()

for (i in 1:length(hla_files)) {
  tmp_file <- hla_files[i]
  print(paste0("Read table: ", tmp_file))
  affinity_table <- read.delim(paste0(OUTPUT, "/", tmp_file), skip = 2, header = F)
  colnames(affinity_table) <- c("Pos", "Peptide_ID", "core", "icore", 
                                "EL_score", "EL_Rank", "BA_score", 
                                "BA_Rank", "Ave", "NB", "N")
  
  cols_to_keep <- c("Peptide_ID", "icore", "BA_score", "BA_Rank")

  tmp_strong <- affinity_table %>%
    filter(BA_Rank != "BA-score") %>%
    mutate(BA_Rank = as.numeric(BA_Rank)) %>%
    filter(BA_Rank <= 0.5) %>%
    select(all_of(cols_to_keep)) %>%
    mutate(
      Sample = gsub("/.*$", "", tmp_file),
      HLA_allele = tmp_file %>%
        gsub("..*HLA", "HLA", .) %>%
        gsub("\\.xls", "", .) %>%
        gsub("_", ":", .)
    )
  all_strong <- rbind(all_strong, tmp_strong)
  
  tmp_weak <- affinity_table %>%
    filter(BA_Rank != "BA-score") %>%
    mutate(BA_Rank = as.numeric(BA_Rank)) %>%
    filter(BA_Rank <= 2 & BA_Rank > 0.5) %>%
    select(all_of(cols_to_keep)) %>%
    mutate(
      Sample = gsub("/.*$", "", tmp_file),
      HLA_allele = tmp_file %>%
        gsub("..*HLA", "HLA", .) %>%
        gsub("\\.xls", "", .) %>%
        gsub("_", ":", .)
    )
  all_weak <- rbind(all_weak, tmp_weak)
}

# Function to calculate TAP efficiency scores for 9-mer peptides
# Consensus matrix downloaded from: https://doi.org/10.4049/jimmunol.171.4.1741
calculate_tap_binding <- function(peptides_df) {
  # Create scoring matrix from the consensus scores
  # Replace question marks with empty strings in the raw data
  scoring_matrix <- matrix(
    c(
      -1.56, -0.25, -0.10,  0.24, -0.10,  0.17,  0.27,  0.00,  0.55,  # A
      0.05, -0.01, -0.02,  0.11,  0.09,  0.05,  0.00, -0.13,  0.00,  # C
      1.37,  1.42,  1.83, -0.23,  0.33,  0.32,  1.07,  0.32,  1.83,  # D
      1.65,  0.02,  1.51,  0.08,  0.54, -0.13,  0.64,  0.44,  1.58,  # E
      1.03, -0.45, -1.05, -0.50, -0.26,  0.08, -0.50,  0.17, -2.52,  # F
      0.28,  1.14,  1.70,  0.45,  0.66,  0.12,  1.41, -0.38,  1.41,  # G
      0.21,  0.33, -0.23, -0.21, -0.11, -0.06, -0.19,  0.39,  0.55,  # H
      -0.11, -0.49, -0.62, -0.09, -0.42, -0.75, -0.94,  0.45, -0.52,  # I
      -1.03, -0.41,  0.09, -0.23, -0.08, -0.26,  0.44,  0.12, -0.45,  # K
      -0.50,  0.09, -0.11,  0.11, -0.34,  0.02, -0.73,  0.01, -0.94,  # L
      -0.38, -0.46, -0.58, -0.35, -0.26,  0.30, -0.64, -0.11, -0.29,  # M
      -1.43,  0.69,  1.01,  0.38,  0.49, -0.27,  0.16,  0.33,  1.33,  # N
      1.43,  3.00,  0.22, -0.04, -0.72, -0.13, -0.84,  0.03, -0.09,  # P
      0.47, -0.97,  0.39,  0.15,  0.15, -0.07,  0.34,  0.26,  0.12,  # Q
      -1.34, -1.47, -0.42, -0.27, -0.32, -0.75, -0.09, -0.42, -1.47,  # R
      -0.56, -0.34,  0.11,  0.27,  0.45,  0.31,  0.87, -0.51,  2.26,  # S
      -0.12, -0.04,  0.43,  0.23,  0.43,  0.49,  0.39, -0.46,  0.72,  # T
      -0.49, -0.50, -0.71,  0.27,  0.37, -0.02, -0.29,  0.10, -0.30,  # V
      0.54, -0.64, -1.65, -0.18, -0.78,  0.31, -0.50, -0.63, -0.87,  # W
      0.50, -0.67, -1.80, -0.18, -0.13,  0.28, -0.87,  0.02, -2.91   # Y
    ),
    nrow = 20,
    byrow = TRUE
  )
  rownames(scoring_matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  colnames(scoring_matrix) <- paste0("pos", 1:9)
  
  # Function to calculate score for a single peptide
  score_peptide <- function(peptide) {
    aa_vector <- strsplit(peptide, "")[[1]]
    if(length(aa_vector) != 9) {
      warning(sprintf("Peptide %s length is not 9", peptide))
      return(NA)
    }
    position_scores <- sapply(1:9, function(i) {
      if(aa_vector[i] %in% rownames(scoring_matrix)) {
        scoring_matrix[aa_vector[i], i]
      } else {0}
    })
    sum(position_scores)
  }
  
  # Calculate TAP efficiency scores for all peptides
  peptides_df$TAP_logIC50 <- sapply(peptides_df$Peptide_ID, score_peptide)
  
  return(peptides_df)
}

write.csv(all_strong, paste0(OUTPUT, "/Strong_binders.csv"))
write.csv(all_weak, paste0(OUTPUT, "/Weak_binders.csv"))

new_strong <- calculate_tap_binding(all_strong) %>%
  arrange(TAP_logIC50, BA_Rank) %>%
  filter(TAP_logIC50 < 0) %>%
  mutate(weight = log2(abs(1/BA_Rank * TAP_logIC50) + 1)) %>%
  select(icore, HLA_allele, weight)
write.csv(new_strong, paste0(OUTPUT, "/Strong_network.csv"))

new_weak <- calculate_tap_binding(all_weak) %>%
  arrange(TAP_logIC50, BA_Rank) %>%
  filter(TAP_logIC50 < 0) %>%
  mutate(weight = log2(abs(1/BA_Rank * TAP_logIC50) + 1)) %>%
  select(icore, HLA_allele, weight)

write.csv(new_weak, paste0(OUTPUT, "/Weak_network.csv"))