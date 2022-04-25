# Generate a table of genes linked to body size and development, ranked by the 
# strength of different types of evidence
#
# Yury V Bukhman
# 21 Dec 2021
# Project: BWGENOME/Gene_lists/1_body_size_and_development

library(tidyverse)

# Directory paths
GENE_LISTS <- "/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/BWGENOME-439_text_mining"
setwd("/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/Gene_lists/1_body_size_and_development")

# Curated gene lists
# - whales
whales <- paste0(GENE_LISTS,"/curated_gene_lists/Lagunas-Rangel_2021_whales.txt") %>% read_lines()
tollis <- paste0(GENE_LISTS,"/curated_gene_lists/Tollis_2019_whales.txt") %>% read_lines()
whales <- c(whales, tollis) %>% unique
gene_table <- tibble(gene = whales, whales_lit = 1)

# - dogs
dogs <- paste0(GENE_LISTS,"/curated_gene_lists/Ostrander_2017_dogs.txt") %>% read_lines() %>% unique
dogs_table <- tibble(gene = dogs, Ostrander_dogs = 1)
gene_table <- merge(gene_table, dogs_table, by = "gene", all.x = TRUE, all.y = TRUE)

# - dwarfism
dwarfs <- paste0(GENE_LISTS,"/curated_gene_lists/gene_disease_db_dwarfism.txt") %>% read_lines() %>% unique
dwarfs_table <- tibble(gene = dwarfs, curated_dwarfism = 1)
gene_table <- merge(gene_table, dwarfs_table, by = "gene", all.x = TRUE, all.y = TRUE)

# - developmental
development <- paste0(GENE_LISTS,"/curated_gene_lists/Developmental_genes_Chu.txt") %>% 
  read_lines %>% str_remove_all("\\s") %>% unique
development_table <- tibble(gene = development, developmental = 1)
gene_table <- merge(gene_table, development_table, by = "gene", all.x = TRUE, all.y = TRUE)

# - cattle
cattle <- paste0(GENE_LISTS,"/curated_gene_lists/Bouwman_2018_cattle.txt") %>% read_lines() %>% unique
cattle_table <- tibble(gene = cattle, Bouwman_cattle = 1)
gene_table <- merge(gene_table, cattle_table, by = "gene", all.x = TRUE, all.y = TRUE)

# - sheep
sheep <- paste0(GENE_LISTS,"/curated_gene_lists/Kominakis_2017_sheep.txt") %>% read_lines() %>% unique
sheep_table <- tibble(gene = sheep, Kominakis_sheep = 1)
gene_table <- merge(gene_table, sheep_table, by = "gene", all.x = TRUE, all.y = TRUE)

# KinderMiner
km_genes <- c() # this vector will contain a union of all genes found by KinderMiner searches

# - Kalpana's files
files <- paste0(GENE_LISTS, "/Kalpana/KM_gene_vs_diseases_BlueWhales/KM_condition_vs_genes/") %>% dir
for (f in files) {
  if (f == "README.txt") next
  km <- paste0(GENE_LISTS, "/Kalpana/KM_gene_vs_diseases_BlueWhales/KM_condition_vs_genes/",f) %>% read_tsv()
  genes <- km$target %>% str_remove(., "\\:.*")
  km_genes <- union(km_genes,genes)
}

# - Ron's files
files <- paste0(GENE_LISTS, "/Ron_2021-12-13/") %>% dir
for (f in files) {
  if (f == "README.txt") next
  km <- paste0(GENE_LISTS, "/Ron_2021-12-13/",f) %>% read_tsv()
  genes <- km$target
  km_genes <- union(km_genes,genes)
}

# - merge KinderMiner genes into gene_table
km_table <- tibble(gene = km_genes, KinderMiner = 1)
gene_table <- merge(gene_table, km_table, by = "gene", all.x = TRUE, all.y = TRUE)

# SKiM
skim_genes <- c() # this vector will contain a union of all genes found by KinderMiner searches
files <- paste0(GENE_LISTS, "/Kalpana/SKiM_for_blue_whales/") %>% dir
for (f in files) {
  if (f == "README.txt") next
  skim <- paste0(GENE_LISTS, "/Kalpana/SKiM_for_blue_whales/",f) %>% read_tsv()
  genes <- skim$target %>% str_remove(., "\\:.*")
  skim_genes <- union(skim_genes,genes)
}

# - merge SKiM genes into gene_table
skim_table <- tibble(gene = skim_genes, SKiM = 1)
gene_table <- merge(gene_table, skim_table, by = "gene", all.x = TRUE, all.y = TRUE)

# Assign scores
weights <- list(whales_lit = 4, Ostrander_dogs = 4, curated_dwarfism = 4, 
                developmental = 3, Bouwman_cattle = 3, Kominakis_sheep = 3,
                KinderMiner = 2, SKiM = 1)
for (colname in names(weights)) {
  gene_table[!is.na(gene_table[[colname]]),colname] <- weights[[colname]]
  gene_table[is.na(gene_table[[colname]]),colname] <- 0
}

# Compute composite score, sort, and save
gene_table$score <- apply(gene_table[names(weights)], 1, sum)
gene_table %>% arrange(desc(score), gene) %>% 
  write_csv("genes_of_interest_scored_by_evidence.csv")
