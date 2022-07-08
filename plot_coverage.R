# Plot PacBio read coverage around genes of interest
# Yury V Bukhman, 07 July 2022

SCAFFOLD <- "NC_045803.1"
GENES <- tibble(gene = c("XRCC1","LOC118885654"),
                start = c(11719062, 11788759),
                end = c(11740692, 11812211))
ZOOM_OUT <- 1 # larger number will zoom out more
AVERAGE_COVERAGE <- 51.16

library(tidyverse)

# Read in the coverage bins
setwd("/Volumes/home/Projects/BWGENOME/final_stage_analyses/SegDups")
coverage_bins <- read_tsv("mBalMus1.pri.coverage.bins", 
                          col_names = c("scaffold","start","end","coverage"))

# Compute genomic interval to plot
interval_midpoint <- (min(GENES$start) + max(GENES$end)) / 2
interval_width <- (max(GENES$end) - min(GENES$start)) * ZOOM_OUT
interval_start <- interval_midpoint - interval_width/2
interval_end <- interval_midpoint + interval_width/2

# Retrieve the coverage data to be plotted
interval_data <- filter(coverage_bins, scaffold == SCAFFOLD, 
                        start > interval_start, end < interval_end)
interval_data <- interval_data %>% mutate(midpoint = (start+end)/2)

# Plot the coverage data
ggplot(data=interval_data) + geom_line(aes(x=midpoint,y=coverage)) + 
  geom_hline(yintercept = AVERAGE_COVERAGE, col = "blue") +
  geom_segment(aes(x=start, xend=end, y=0.2*AVERAGE_COVERAGE, yend=0.2*AVERAGE_COVERAGE),
               data = GENES, col="red", lwd=2) +
  geom_text(aes(x=start, y=0.2*AVERAGE_COVERAGE, label=gene), data = GENES, 
            col="red", hjust="left", vjust=1.5)

