#!/usr/bin/env Rscript

# Using a multiple alignment, find sites within a gene that correlate with body
# size. Estimate FDR using a permutation test

# 2-factor model accounting for the habitat

# Yury V Bukhman, 01 Mar 2022 
# Project: /BWGENOME/Interesting_genes/IGF1/body-size-analysis/
#          body-size-loci_try1_cetaceans-only/4_permut-fdr

# Function to compute p values =================================================
compute_pvals <- function(alignment_matrix, my_species, body_mass, bw_start_col, 
                          bw_stop_col, all_start_col, all_stop_col) {
  pvals <- double()
  for (i in bw_start_col:bw_stop_col) {
    # Report progress
    if (i %% 1e4 == 0) cat(i, "\n")
    
    # Build a data frame of bases and weights
    bases <- data.frame(species = my_species, base = alignment_matrix[,i])
    data1 <- merge(body_mass, bases, by = "species")
    
    # If we are near the edges of the alignment matrix, filter out missing species
    if (i < all_start_col | i > all_stop_col) data1 <- filter(data1, base != "-")
    
    # Compute the weight vs. base p value
    base_counts <- summary(as.factor(data1$base))
    if (nrow(data1) == 0 | 
        length(unique(data1$habitat)) == 1 | 
        length(base_counts) == 1 | 
        sort(base_counts, decreasing = TRUE)[2] < 2) {
      # Can't do statistics unless all of the following conditions are met:
      # 1. We have at least some non-missing data
      # 2. Both habitats are represented
      # 3. We have at least 2 kinds of bases
      # 4. We have at least 2 instances of the minor allele
      pval <- 1
    } else {
      # Compute p value using a linear model
      lm1 <- lm(log(adult_body_mass) ~ habitat + base, data = data1)
      lmsum <- summary(lm1)
      if (nrow(lmsum$coefficients) == 3) {
        pval <- lmsum$coefficients[3,"Pr(>|t|)"]
      } else {
        # the base is confounded with the habitat
        pval <- 1
      }
    }
    
    # Append p value to the vector of p values
    pvals <- c(pvals, pval)
  }
  return(pvals)
}

# MAIN =========================================================================
# Libraries
cat("Loading libraries...\n")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sequencing))
suppressPackageStartupMessages(library(optparse))

# Specify command line options
option_list <- list( 
  make_option(c("-a", "--alignment"), help = "multiple alignment file"),
  make_option(c("-p", "--pantheria"), help = "Pantheria database file [default %default]",
              default = "/w5home/ybukhman/Projects/BWGENOME/BWGENOME-337_bodysize_vs_genelength/Data_analysis/PanTHERIA_1-0_WR05_Aug2008.txt"),
  make_option(c("-n", "--num_permutations"), default = 100, 
              help = "number of permutations to estimate FDR [default %default]")
)

# Get command line options 
opt <- parse_args(OptionParser(option_list=option_list))

# Read in the multiple sequence alignment
origMAlign <- readDNAMultipleAlignment(filepath = opt$alignment, format="clustal")

# Extract the alignment matrix
alignment_matrix <- as.matrix(origMAlign@unmasked)

# Properties of the alignment matrix ==========================================
cat("\nProperties of the alignment matrix:\n")

# Define the range of the alignment columns that are present in all species. It 
# starts at the first position that has no "-" symbols and stops at the last 
# such position

# Scan forward from the beginning of the matrix to find the start
for (i in 1:ncol(alignment_matrix)) {
  if (sum(alignment_matrix[,i] == "-") == 0) break
}
all_start_col <- i

# Scan backward from the end of the matrix to find the stop
for (i in ncol(alignment_matrix):1) {
  if (sum(alignment_matrix[,i] == "-") == 0) break
}
all_stop_col <- i

cat("All species are present between columns", all_start_col, "and", 
    all_stop_col, "\n")

# Define the range of the alignment columns that are present in blue whale
# Scan forward from the beginning of the matrix to find the start
for (i in 1:ncol(alignment_matrix)) {
  if (alignment_matrix["Balaenoptera_musculus",i] != "-") break
}
bw_start_col <- i

# Scan backward from the end of the matrix to find the stop
for (i in ncol(alignment_matrix):1) {
  if (alignment_matrix["Balaenoptera_musculus",i] != "-") break
}
bw_stop_col <- i

cat("blue whale range is", bw_start_col, "-", bw_stop_col, "\n")

# Body mass values =============================================================
cat("\nSetting body mass values...\n")

# Read in PanTHERIA database
# col_types = cols() suppresses the column specification message
pantheria <- read_tsv(opt$pantheria, col_types = cols())

# Set adult body mass values
my_species <- rownames(alignment_matrix) %>% str_replace_all("_", " ")
body_mass <- pantheria[pantheria$MSW05_Binomial %in% my_species, c("MSW05_Binomial","5-1_AdultBodyMass_g")]
names(body_mass) <- c("species","adult_body_mass")

# Manually add species missing in Pantheria
body_mass$adult_body_mass[body_mass$species == "Bos taurus"] <- 816466
body_mass <- rbind(body_mass, data.frame(species = "Bos indicus", adult_body_mass = 1e6))
body_mass <- rbind(body_mass, data.frame(species = "Bos indicus x Bos taurus", adult_body_mass = 9e5))
body_mass <- rbind(body_mass, data.frame(species = "Bos mutus", adult_body_mass = 1.2e6))
body_mass <- rbind(body_mass, data.frame(species = "Bison bison bison", adult_body_mass = subset(pantheria, MSW05_Binomial == "Bison bison")$`5-1_AdultBodyMass_g`))
body_mass <- rbind(body_mass, data.frame(species = "Cervus canadensis", adult_body_mass = 497e3))
body_mass <- rbind(body_mass, data.frame(species = "Camelus ferus", adult_body_mass = subset(pantheria, MSW05_Binomial == "Camelus bactrianus")$`5-1_AdultBodyMass_g`))
body_mass <- rbind(body_mass, data.frame(species = "Neophocaena asiaeorientalis as", adult_body_mass = 71800))
body_mass <- rbind(body_mass, data.frame(species = "Balaenoptera acutorostrata sca", adult_body_mass = subset(pantheria, MSW05_Binomial == "Balaenoptera acutorostrata")$`5-1_AdultBodyMass_g`))
body_mass <- rbind(body_mass, data.frame(species = "Vicugna pacos", adult_body_mass = 65e3))
body_mass <- rbind(body_mass, data.frame(species = "Odocoileus virginianus texanus", adult_body_mass = subset(pantheria, MSW05_Binomial == "Odocoileus virginianus")$`5-1_AdultBodyMass_g`))

# Drop the human
body_mass <- body_mass %>% filter(species != "Homo sapiens")

# Assign habitat
body_mass <- body_mass %>% mutate(habitat = "terrestrial")
body_mass$habitat[body_mass$species %in% c("Balaenoptera musculus",
                                           "Delphinapterus leucas",
                                           "Globicephala melas",
                                           "Lagenorhynchus obliquidens",
                                           "Monodon monoceros",
                                           "Lipotes vexillifer",
                                           "Orcinus orca",
                                           "Phocoena sinus",
                                           "Physeter catodon",
                                           "Neophocaena asiaeorientalis as",
                                           "Balaenoptera acutorostrata sca")] <- "aquatic"

# Linear models on the columns of the alignment matrix =========================
cat("\nComputing p values...\n")
pvals <- compute_pvals(alignment_matrix, my_species, body_mass, bw_start_col, 
                       bw_stop_col, all_start_col, all_stop_col)

# Compute how many sites pass various p value thresholds
num_hits <- data.frame(threshold = seq(from=1, to=6, by=0.1), num=NA)
num_hits$num <- apply(num_hits, 1, function(x) sum(-log10(pvals) > x["threshold"], na.rm = TRUE))

# Resampling====================================================================
cat("\nResampling...\n")
resamp_num_hits <- data.frame(threshold = seq(from=1, to=6, by=0.1))
for (iteration in 1:opt$num_permutations) {
  cat("iteration",iteration,"\n")
  
  # Resample body mass values
  resamp_body_mass <- body_mass
  terr <- resamp_body_mass$habitat == "terrestrial"
  resamp_body_mass$adult_body_mass[terr] <- sample(body_mass$adult_body_mass)[terr]
  resamp_body_mass$adult_body_mass[!terr] <- sample(body_mass$adult_body_mass)[!terr]
  
  # Compute p values
  iter_pvals <- compute_pvals(alignment_matrix, my_species, resamp_body_mass, bw_start_col, 
                              bw_stop_col, all_start_col, all_stop_col)
  
  # Compute how many sites pass various p value thresholds
  iter_num_hits <- apply(resamp_num_hits, 1, function(x) sum(-log10(iter_pvals) > x["threshold"]))
  resamp_num_hits <- cbind(resamp_num_hits, iter_num_hits)
}

# FDR and quantile estimation===================================================
cat("\nEstimating FDR and quantile...\n")
# Compute FDR
num_hits$resamp_mean_hits <- apply(resamp_num_hits[,-1], 1, mean)
num_hits$fdr_est <- num_hits$resamp_mean_hits / num_hits$num

# Compute quantile
num_hits$quantile <- NA
for (i in 1:nrow(num_hits)) {
  num_hits$quantile[i] <- sum(resamp_num_hits[i,-1] < num_hits$num[i]) / 
    (ncol(resamp_num_hits) - 1)
}

# Write out results=============================================================
cat("\nWriting out results...\n")
write_csv(tibble(alignment_column = bw_start_col:bw_stop_col, pval = pvals), 
          "p_values.csv")
write_csv(num_hits, "num_hits_and_fdr.csv")
save.image("R-workspace-image.RData")