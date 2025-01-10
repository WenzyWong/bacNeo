#################################################################################
# Rscript for extracting strong and weak binders from netMHCpan output files

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-01-10
#################################################################################
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
OUTPUT <- args[1]

hla_files <- list.files(path = OUTPUT, pattern = "HLA-", recursive = T)

all_strong <- data.frame()
all_weak <- data.frame()

for (i in 1:length(hla_files)) {
  tmp_file <- hla_files[i]
  affinity_table <- read.delim(paste0(OUTPUT, "/", file), skip = 2, header = F)
  colnames(affinity_table) <- c("Pos", "Peptide_ID", "core", "icore", 
                                "EL_score", "EL_Rank", "BA_score", 
                                "BA_Rank", "Ave", "NB", "N")
  
  cols_to_keep <- c("Peptide_ID", "icore", "BA_score", "BA_Rank")
  tmp_strong <- affinity_table %>%
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
write.csv(all_strong, paste0(OUTPUT, "/Strong_binders.csv"))
write.csv(all_weak, paste0(OUTPUT, "/Weak_binders.csv"))