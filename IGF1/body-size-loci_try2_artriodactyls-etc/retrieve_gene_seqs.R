# Retrieve genomic sequences of a gene + flanking regions for all members of a taxon
# Yury V Bukhman, 2022-02-04
# Project: BWGENOME/Interesting_genes/IGF1

# Initialize
library(rentrez)
setwd("/Users/ybukhman/Google Drive/Morgridge/Projects/BWGENOME/Interesting_genes/IGF1/body-size-loci_try2_artriodactyls-etc/1_genomic_seqs")
GENE <- "IGF1"
TAXA <- c("Artiodactyla","Canis lupus familiaris","Homo sapiens")
FLANK <- 1e4 # the length of the flanking regions arround the gene to retrieve

for (TAXON in TAXA) {
  # Get all instances of a specified gene in a taxon from NCBI Gene database
  term <- paste0(TAXON, "[Organism] AND ", GENE, "[Gene Name]")
  genes <- entrez_search(db = "gene", term = term, retmax = 40)
  genes_sum <- entrez_summary(db = "gene", id = genes$ids)
  if (length(genes$ids) == 1) genes_sum <- list(genes_sum)
  
  # Retrieve sequences and write to a file
  output_file <- paste0(GENE, "_genomic_seqs.fna")
  for (gene in (genes_sum)) {
    if (gene$organism$scientificname == "Tursiops truncatus" | length(gene$locationhist) == 0) next
    cat(gene$organism$scientificname, "\n")
    
    # Get coordinates of the sequence to retrieve
    seq_accession <- gene$locationhist$chraccver[1]
    gene_start <- gene$locationhist$chrstart[1]
    gene_stop <- gene$locationhist$chrstop[1]
    if (gene_start < gene_stop) {
      strand <- 1
      seq_start <- gene_start - FLANK
      seq_stop <- gene_stop + FLANK
    } else {
      strand <- 2
      seq_stop <- gene_start + FLANK
      seq_start <- gene_stop - FLANK
    }
    
    # Fetch the sequence
    seq <- entrez_fetch(db = "nuccore", id = seq_accession, rettype="fasta", 
                        strand = strand, seq_start = seq_start, seq_stop = seq_stop)
    
    # Replace the header with the scientific name of the organism
    org <- str_replace_all(gene$organism$scientificname, " ", "_")
    seq <- str_replace(seq, "^>.+\n", paste0(">",org,"\n"))
    
    # Write out to a file
    write_file(seq, output_file, append = TRUE)
  }
}