#################################################################################
# Rscript for extracting peptides with species information
# and draw relation-network

# Yunzhe WANG, yunzhewang24@m.fudan.edu.cn
# Updated: 2025-01-21
#################################################################################
if(!require("Rcpp")) {
    install.packages("Rcpp", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}
if(!require("withr")) {
  install.packages("withr", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}
if (!require("readxl")) {
  withr::with_makevars(c(PKG_LIBS = "-liconv"), install.packages("readxl", repo = "https://mirrors.ustc.edu.cn/CRAN/"), assignment = "+=")
}
if(!require("stringr")) {
  install.packages("stringr", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}
if(!require("data.table")) {
  install.packages("data.table", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}
if(!require("dplyr")) {
  install.packages("dplyr", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}
if(!require("igraph")) {
  install.packages("igraph", repo = "https://mirrors.ustc.edu.cn/CRAN/")
}

library(readxl)
library(stringr)
library(igraph)

args = commandArgs(trailingOnly=TRUE)
OUTPUT <- args[1]

find_bac <- function (description, infoType) {
  library(stringr)
  library(dplyr)
  genusInfo <- c()
  speciesInfo <- c()
  for (i in c(1 : length(description))) {
    tmp <- unlist(strsplit(description[i], ";"))
    len <- length(tmp)
    genus_tmp <- c()
    species_tmp <- c()
    for (j in 1:len) {
      genus_tmp <- c(genus_tmp, 
                    unlist(str_split(unlist(strsplit(tmp[j], "="))[2], " ", 3))[1])
      species_tmp <- c(species_tmp, 
                      paste(unlist(str_split(unlist(strsplit(tmp[j], "="))[2], " ", 3))[1],
                            unlist(str_split(unlist(strsplit(tmp[j], "="))[2], " ", 3))[2]))
    }
    genus_tmp <- unique(genus_tmp)
    if (length(genus_tmp) > 1) {genus_tmp <- paste(genus_tmp, collapse = ";")}
    genusInfo <- c(genusInfo, genus_tmp)
    species_tmp <- unique(species_tmp)
    if (length(species_tmp) > 1) {species_tmp <- paste(species_tmp, collapse = ";")}
    speciesInfo <- c(speciesInfo, species_tmp)

    res <- case_when(
      infoType == "g" ~ genusInfo,
      infoType == "s" ~ speciesInfo
    )
  }
  return(res)
}

FILE <- list.files(path = OUTPUT, pattern = "proteinGroups.txt",
                   recursive = T)
print(gsub("/combined/txt/proteinGroups.txt", "", FILE[1]))
tmpDf <- data.table::fread(paste0(OUTPUT, "/", FILE[1]),
                           header = T, data.table = F,
                           check.names = F)
# Removing contaminants, reversive peptides and human-derived peptides
tmpDf <- tmpDf[tmpDf$`Potential contaminant` != "+" & tmpDf$Reverse != "+", ]
tmpDf <- tmpDf[!grepl("_HUMAN", tmpDf$`Fasta headers`), ]
tmpDf <- tmpDf[ , grepl("Unique peptides |Protein IDs|Fasta headers|Peptide sequences", 
                        colnames(tmpDf))]

# Combining different names which were assigned to the same peptides
# Using ";" as combining mark and adjusting the names in alphabetic order
tmpName <- strsplit(tmpDf$`Protein IDs`, ";")
reName <- c()
for (i in 1:length(tmpName)) {
  reName <- c(reName, paste(sort(tmpName[[i]]), collapse = ";"))
}
tmpDf <- cbind(Protein.ID = reName, tmpDf[ , -1])
Df <- tmpDf
print(paste("Unique peptide numbers:", nrow(tmpDf)))

for (f in FILE[2:length(FILE)]) {
  print(gsub("/combined/txt/proteinGroups.txt", "", f))
  
  tmpDf <- data.table::fread(paste0(OUTPUT, "/", f),
                             header = T, data.table = F,
                             check.names = F)
  tmpDf <- tmpDf[tmpDf$`Potential contaminant` != "+" & 
                   tmpDf$Reverse != "+", ]
  tmpDf <- tmpDf[!grepl("_HUMAN", tmpDf$`Fasta headers`), ]
  tmpDf <- tmpDf[ , grepl("Unique peptides |Protein IDs|Fasta headers|Peptide sequences", 
                          colnames(tmpDf))]
  
  tmpName <- strsplit(tmpDf$`Protein IDs`, ";")
  reName <- c()
  for (i in 1:length(tmpName)) {
    reName <- c(reName, paste(sort(tmpName[[i]]), collapse = ";"))
  }
  tmpDf <- cbind(Protein.ID = reName, tmpDf[ , -1])
  print(paste("Unique peptide numbers:", nrow(tmpDf)))
  
  Df <- merge(Df, tmpDf, by = c("Protein.ID", "Fasta headers", "Peptide sequences"), all = T)
}
Df[is.na(Df)] <- 0
# dim(Df)
aggDf <- aggregate(Df[ , 4:ncol(Df)], by = list(Df$Protein.ID), FUN = sum)
aggDf <- merge(Df[1:3], aggDf, by.x = "Protein.ID", by.y = "Group.1", all.y = T)
for (i in unique(aggDf$Protein.ID[duplicated(aggDf$Protein.ID)])) {
  seqs <- aggDf$`Peptide sequences`[aggDf$Protein.ID == i]
  seqs <- seqs[nchar(seqs) == max(nchar(seqs))][1]
  aggDf$`Peptide sequences`[aggDf$Protein.ID == i] <- seqs
}
aggDf <- aggDf[!duplicated(aggDf$Protein.ID), ]
# dim(aggDf)
aggDf$Genus <- find_bac(aggDf$`Fasta headers`, "g")
aggDf$Species <- find_bac(aggDf$`Fasta headers`, "s")
write.csv(aggDf, paste0(OUTPUT, "/peptide_info.csv"))
seq <- unlist(strsplit(aggDf$`Peptide sequences`, ";"))
seq <- seq[nchar(seq) >= 8]
write.table(seq, row.names = F, paste0(OUTPUT, "/sequences.txt"))

allsp <- unlist(strsplit(aggDf$Species, ";"))
allpro <- c()
allline <- c()
for (i in 1:nrow(aggDf)) {
  if (grepl(";", aggDf$Species[i])) {
    reptime <- length(unlist(strsplit(aggDf$Species[i], ";")))
    allpro <- c(allpro, rep(aggDf$Protein.ID[i], reptime))
    allline <- c(allline, rep(sum(aggDf[i, 4:ncol(Df)]), reptime))
  } else {
    allpro <- c(allpro, aggDf$Protein.ID[i])
    allline <- c(allline, sum(aggDf[i, 4:ncol(Df)]))
  }
}
link <- data.frame(
  node1 = allpro,
  node2 = allsp,
  linesize = allline
)

c(as.character(link$node1), as.character(link$node2)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")
head(vertices)

g <- graph_from_data_frame(link, vertices = vertices, directed = FALSE)

E(g)$color <- "grey30"
E(g)$width <- E(g)$linesize

V(g)$color <- if_else(grepl(" ", V(g)$name), "burlywood", "darkseagreen")
V(g)$size <- if_else(grepl(" ", V(g)$name), 5, 3)
V(g)$label.cex <- if_else(grepl(" ", V(g)$name), 0.5, 0.3)
V(g)$label.color <- if_else(grepl(" ", V(g)$name), "orange4", "darkolivegreen")

top_edges <- order(E(g)$linesize, decreasing = TRUE)[1:3]

top_edge_groups <- lapply(top_edges, function(e) {
  ends <- ends(g, e)  # Get the two vertices of the edge
  c(ends[1], ends[2]) # Return vertex indices
})

pdf(file.path(OUTPUT, "/species_peptide_network.pdf"), width = 10, height = 10)
print(plot.igraph(g, 
                  edge.curved = 0,
                  vertex.frame.color = "grey",
                  vertex.label.dist = 0,
                  mark.groups = top_edge_groups,
                  mark.col = NA,
                  mark.border = "red",
                  mark.lty = 2
))
dev.off()
