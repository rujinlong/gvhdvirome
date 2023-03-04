#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)
library(phyloseq)

# Add command line arguments
p <- arg_parser("Phyloseq analysis of virome data")
p <- add_argument(p, "--metadata", help="metadata table")
p <- add_argument(p, "--votu", help="vOTU table")
p <- add_argument(p, "--taxa", help="vTaxa table")
p <- add_argument(p, "--output", help="Output phyloseq object (.RData)")
argv <- parse_args(p)

vmeta <- read_tsv(argv$metadata) %>% 
    column_to_rownames("sample_id_virome") %>% 
    sample_data()

votu <- read_delim(argv$votu, delim = ",") %>% 
    column_to_rownames("contig_id") %>% 
    as.matrix() %>% 
    otu_table(taxa_are_rows = TRUE)

vtax <- read_delim(argv$taxa, delim = ",") %>% 
    column_to_rownames("contig_id") %>% 
    as.matrix() %>% 
    tax_table()

# Create phyloseq object
pseq_vir <- phyloseq(votu, vtax, vmeta)

save(pseq_vir, file=argv$output)
