# Generate a table of genes linked to body size, ranked by the 
# strength of different types of evidence
#
# Yury V Bukhman
# 21 Jan 2021
# Project: BWGENOME/Gene_lists/4_body_size

library(tidyverse)

# Directory paths
GENE_LISTS <- "/Volumes/home/Projects/BWGENOME/BWGENOME-439_text_mining"
DIFF_COPY_NUMBER_GENES <- "/Volumes/home/Projects/BWGENOME/final_stage_analyses/Gene_lists/2_differential_copy_number/Version_2/diff_copy_num_genes.csv"
setwd("/Volumes/home/Projects/BWGENOME/final_stage_analyses/Gene_lists/4_body_size/Version_2")

# Curated gene lists
# - dogs
dogs <- paste0(GENE_LISTS,"/curated_gene_lists/Ostrander_2017_dogs.txt") %>% read_lines() %>% unique
gene_table <- tibble(gene = dogs, Ostrander_dogs = 1)

# - dwarfism
dwarfs <- paste0(GENE_LISTS,"/curated_gene_lists/gene_disease_db_dwarfism.txt") %>% read_lines() %>% unique
dwarfs_table <- tibble(gene = dwarfs, curated_dwarfism = 1)
gene_table <- merge(gene_table, dwarfs_table, by = "gene", all.x = TRUE, all.y = TRUE)

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
files <- c("Body_size/C0005901_body_size_evaluationResult_defaultpvaue.txt",
           "Dwarfism/C0013336_dwarfism_evaluationResult_defaultpvaue.txt",
           "Gigantism/C0017547_Genetic_giant_evaluationResult_defaultpvaue.txt",
           "Overgrowth/C1849265_overgrowth_evaluationResult_defaultpvalue.txt")
for (f in files) {
  if (f == "README.txt") next
  km <- paste0(GENE_LISTS, "/Kalpana/KM_gene_vs_diseases_BlueWhales/KM_condition_vs_genes/",f) %>% read_tsv()
  genes <- km$target %>% str_remove(., "\\:.*")
  km_genes <- union(km_genes,genes)
}

# - Ron's file
file <- "results_control_body_size_kinderminer_results_match_all_pval_lt_1e-5 - 2021-12-13T152714.287.tsv"
km <- paste0(GENE_LISTS, "/Ron_2021-12-13/",file) %>% read_tsv()
genes <- km$target
km_genes <- union(km_genes,genes)

# - merge KinderMiner genes into gene_table
km_table <- tibble(gene = km_genes, KinderMiner = 1)
gene_table <- merge(gene_table, km_table, by = "gene", all.x = TRUE, all.y = TRUE)

# Not doing SKiM: too many genes in the SKiM lists

# Assign scores
weights <- list(Ostrander_dogs = 3, curated_dwarfism = 3, 
                Bouwman_cattle = 2, Kominakis_sheep = 2,
                KinderMiner = 1)
for (colname in names(weights)) {
  gene_table[!is.na(gene_table[[colname]]),colname] <- weights[[colname]]
  gene_table[is.na(gene_table[[colname]]),colname] <- 0
}

# Compute composite score, sort, and save
gene_table$score <- apply(gene_table[names(weights)], 1, sum)
gene_table %>% arrange(desc(score), gene) %>% 
  write_csv("body_size_genes_scored_by_evidence.csv")

# Join to differentially duplicated genes
diff_dups <- read_csv(DIFF_COPY_NUMBER_GENES)
ddbs <- merge(diff_dups, gene_table, by = "gene")
ddbs %>% arrange(desc(Blue_whale), desc(score), gene) %>% 
  write_csv("differentially_duplicated_body_size_genes.csv")
