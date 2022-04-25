# Find genes that are lost in blue whale but still present in vaquita or vice
# versa. Intersect this list with genes of potential interest linked 
# to developmental clock, body size, or longevity
#
# Yury V Bukhman
# 29 Dec 2021
# Project: BWGENOME/Gene_lists/3_differential_gene_loss

library(tidyverse)

# Parameters
WORK_DIR <- "/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/Gene_lists/3_differential_gene_loss"
TOGA_ROOT_DIR <- "/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/TOGA/data"
GENES_OF_INTEREST_FILE <- "/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/Gene_lists/1_body_size_and_development/genes_of_interest_scored_by_evidence.csv"
setwd(WORK_DIR)

# Find the genes that are lost in each assembly
find_missing_genes <- function(input_file) {
  orthologs <- read_tsv(paste0(TOGA_ROOT_DIR,input_file))
  orthologs_miss <- filter(orthologs, orthology_class == "one2zero")
  orthologs_miss$t_transcript %>% str_replace("^ENST\\d+\\.","") %>% unique
}
bw_primary <- find_missing_genes("/human_hg38_blue_whale/orthologsClassification.tsv")
bw_alt <- find_missing_genes("/human_hg38_blue_whale_alt/orthologsClassification.tsv")
vaq_primary <- find_missing_genes("/human_hg38_vaquita/orthologsClassification.tsv")
vaq_alt <- find_missing_genes("/human_hg38_vaquita_alt/orthologsClassification.tsv")

# Find the genes that are lost in one species but not the other
bw_homozyg_miss <- intersect(bw_primary, bw_alt)
vaq_homozyg_miss <- intersect(vaq_primary, vaq_alt)
missing_in_bw_not_vaq <- setdiff(bw_homozyg_miss, vaq_homozyg_miss)
missing_in_vaq_not_bw <- setdiff(vaq_homozyg_miss, bw_homozyg_miss)

# - save the gene lists
write_lines(bw_homozyg_miss, "bw_homozyg_miss.txt")
write_lines(vaq_homozyg_miss, "vaq_homozyg_miss.txt")
write_lines(missing_in_bw_not_vaq, "missing_in_bw_not_vaq.txt")
write_lines(missing_in_vaq_not_bw, "missing_in_vaq_not_bw.txt")

# Genes of potential interest missing in one species but not the other
genes_of_interest <- read_csv(GENES_OF_INTEREST_FILE)
genes_of_interest %>% filter(gene %in% missing_in_bw_not_vaq) %>%
  write_csv("genes_of_interest_missing_in_bw_not_vaq.csv")
genes_of_interest %>% filter(gene %in% missing_in_vaq_not_bw) %>%
  write_csv("genes_of_interest_missing_in_vaq_not_bw.csv")

