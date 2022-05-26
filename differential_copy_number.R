# Find genes whose copy numbers differ by a specified fold change between blue 
# whale and vaquita. Intersect this list with genes of potential interest linked 
# to developmental clock, body size, or longevity
#
# Yury V Bukhman
# 29 Dec 2021
# Project: BWGENOME/Gene_lists/2_differential_copy_number

library(tidyverse)

# Parameters
WORK_DIR <- "/Volumes/home/Projects/BWGENOME/final_stage_analyses/Gene_lists/Version_2/2_differential_copy_number"
GENE_COPY_NUMBERS_FILE <- "GeneDups.tsv"
GENES_OF_INTEREST_FILE <- "/Volumes/home/Projects/BWGENOME/final_stage_analyses/Gene_lists/1_body_size_and_development/genes_of_interest_scored_by_evidence.csv"
FOLD_CHANGE_CUTOFF <- 2
setwd(WORK_DIR)

# Find genes that satisfy the fold change cutoff
segdups <- read_tsv(GENE_COPY_NUMBERS_FILE)
fold_change <- (segdups$Blue_whale+2)/(segdups$Vaquita+2)
diff_copy_num <- abs(log2(fold_change)) >= abs(log2(FOLD_CHANGE_CUTOFF))
diffdups <- segdups[diff_copy_num,]
write_csv(diffdups, "diff_copy_num_genes.csv")

# Intersect with genes of interest
genes_of_interest <- read_csv(GENES_OF_INTEREST_FILE)
diff_dup_genes_of_interest <- merge(diffdups, genes_of_interest, by = "gene")
diff_dup_genes_of_interest %>% arrange(desc(score)) %>% write_csv("diff_copy_num_genes_of_interest.csv")
