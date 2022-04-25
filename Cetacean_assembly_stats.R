# Retrieve and plot cetacean assembly stats from NCBI
# Yury V Bukhman
# Projects/BWGENOME/final_stage_analyses/Cetacean_assembly_stats/
# 30 Mar 2022

library(rentrez)
library(tidyverse)
library(ggplot2)

# Assembly records
# Use reference and representative assemblies only. Otherwise, I get 42 assemblies
# that represent only 16 species. 20 of the 42 are different versions of house mouse
r_search <- entrez_search(db="assembly", retmax=200,
                          term='"Cetacea"[Organism] AND (latest[filter] AND ("reference genome"[filter] OR "representative genome"[filter]) AND (all[filter] NOT anomalous[filter] AND all[filter] NOT partial[filter]))')

# Document summaries
multi_summs <- entrez_summary(db="assembly", id=r_search$ids)

# Extract contiguity stats
contiguity_stats <- extract_from_esummary(multi_summs, 
                                          c("speciesname",
                                            "assemblyname",
                                            "assemblyaccession",
                                            "contign50",
                                            "scaffoldn50"))
contiguity_stats <- t(contiguity_stats)  %>% apply(., 2, unlist) %>% as_tibble()
contiguity_stats <- mutate(contiguity_stats, 
                           contign50 = as.integer(contign50),
                           scaffoldn50 = as.integer(scaffoldn50))

# Plot contiguity stats with ContigN50 as X
ggplot(contiguity_stats) + 
  geom_point(aes(x = contign50, y = scaffoldn50)) +
  scale_x_log10() + scale_y_log10() +
  geom_point(data = filter(contiguity_stats, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
            aes(y = scaffoldn50, x = contign50), color = "blue") +
  geom_text(data = filter(contiguity_stats, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
            aes(y = scaffoldn50, x = contign50, label = speciesname, angle = 90), hjust = 1.1, size = 2, color = "blue") +
  geom_point(data = filter(contiguity_stats, speciesname == "Balaenoptera musculus"),
             aes(y = scaffoldn50, x = contign50), color = "red") +
  geom_text(data = filter(contiguity_stats, speciesname == "Balaenoptera musculus"),
            aes(y = scaffoldn50, x = contign50, label = speciesname, angle = 90), hjust = 1.1, size = 2, color = "red") +
  labs(x = "Contig N50", y = "Scaffold N50") +
  theme(text = element_text(family = "Helvetica", size = 8))

#- save the plot
ggsave("Contiguity_stats_3in.pdf", width = 3, height = 3*0.618)

# Extract BUSCO scores
busco_colnames <- multi_summs[[1]]$busco  %>% names
busco_scores <- matrix(nrow = length(multi_summs), ncol = length(busco_colnames), 
                       dimnames = list(extract_from_esummary(multi_summs,"speciesname"),
                                       busco_colnames))
for(i in 1:length(multi_summs)){
  busco_scores[i,] <- multi_summs[[i]]$busco  %>% unlist
}

# Transform to a tibble and toss out empty rows
busco_scores <- as_tibble(busco_scores, rownames = "speciesname")
busco_scores <- filter(busco_scores, complete != "")
busco_scores <- mutate(busco_scores, 
                       complete = as.numeric(complete)*100, 
                       singlecopy = as.numeric(singlecopy)*100,
                       duplicated = as.numeric(duplicated)*100,
                       fragmented = as.numeric(fragmented)*100,
                       missing = as.numeric(missing)*100,
                       totalcount = as.numeric(totalcount))

# Plot BUSCO scores, Complete and Duplicated
ggplot(busco_scores) + 
  geom_point(aes(x = complete, y = duplicated)) +
  geom_point(data = filter(busco_scores, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
            aes(x = complete, y = duplicated), color = "blue") +
  geom_text(data = filter(busco_scores, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
            aes(x = complete, y = duplicated, label = speciesname), hjust = 1.1, size = 2, color = "blue") +
  geom_point(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
             aes(x = complete, y = duplicated), color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
            aes(x = complete, y = duplicated, label = speciesname), hjust = 1.1, size = 2, color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Physeter catodon"),
            aes(x = complete, y = duplicated, label = speciesname), hjust = -0.1, size = 2) +
  labs(x = "% Complete", y = "% Duplicated") +
  theme(text = element_text(family = "Helvetica", size = 8))

#- save the plot
ggsave("BUSCO_scores_CD_3in.pdf", width = 3, height = 3*0.618)

# Plot BUSCO scores, Fragmented and Missing
ggplot(busco_scores) + 
  geom_point(aes(x = fragmented, y = missing)) +
  geom_point(data = filter(busco_scores, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
             aes(x = fragmented, y = missing), color = "blue") +
  geom_text(data = filter(busco_scores, speciesname == "Phocoena sinus"),
            aes(x = fragmented, y = missing, label = speciesname), hjust = -0.1, 
            size = 2, color = "blue") +
  geom_text(data = filter(busco_scores, speciesname == "Tursiops truncatus"),
            aes(x = fragmented, y = missing, label = speciesname), hjust = -0.1, 
            vjust = -0.5, size = 2, color = "blue") +
  geom_point(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
             aes(x = fragmented, y = missing), color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
            aes(x = fragmented, y = missing, label = speciesname), hjust = 0, 
            vjust = 1.7, size = 2, color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Physeter catodon"),
            aes(x = fragmented, y = missing, label = speciesname), hjust = 1.1, size = 2) +
  labs(x = "% Fragmented", y = "% Missing") +
  theme(text = element_text(family = "Helvetica", size = 8))

#- save the plot
ggsave("BUSCO_scores_FM_3in.pdf", width = 3, height = 3*0.618)

# Plot BUSCO scores, Complete and Fragmented
ggplot(busco_scores) + 
  geom_point(aes(x = complete, y = fragmented)) +
  geom_point(data = filter(busco_scores, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
             aes(x = complete, y = fragmented), color = "blue") +
#  geom_text(data = filter(busco_scores, speciesname %in% c("Phocoena sinus", "Tursiops truncatus")),
#            aes(x = complete, y = fragmented, label = speciesname), hjust = 1.1, size = 2, color = "blue") +
  geom_point(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
             aes(x = complete, y = fragmented), color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Balaenoptera musculus"),
            aes(x = complete, y = fragmented, label = speciesname, angle = 90), 
            hjust = -0.1, size = 2, color = "red") +
  geom_text(data = filter(busco_scores, speciesname == "Delphinapterus leucas"),
            aes(x = complete, y = fragmented, label = speciesname, angle = 90), 
            hjust = -0.1, size = 2) +
  geom_text(data = filter(busco_scores, speciesname == "Physeter catodon"),
            aes(x = complete, y = fragmented, label = speciesname), hjust = -0.1, size = 2) +
  labs(x = "% Complete", y = "% Fragmented") +
  theme(text = element_text(family = "Helvetica", size = 8))

#- save the plot
ggsave("BUSCO_scores_CF_3in.pdf", width = 3, height = 3*0.618)
