library(tidyverse)
library(readxl)
library(phyloseq)
library(here)
source(here("scripts/scripts.R"))

# --------- File path --------------
fpath_metadata = here("data/processed/metadata.tsv")
fpath_16s <- here("data/collab/2021-10-20_MOFA-additional-notes_email/TaxaSubset_Genus_V1V3.xlsx")
fpath_ITS <- here("data/collab/2021-10-20_MOFA-additional-notes_email/TaxaSubset_Genus_ITS.xlsx")
# fpath_vabundance = here("data/hms/virmerge/files/abundance/abundance_contigs_count.tsv")
fpath_vabundance = here("vabundance.tsv")
fpath_vOTUs = here("data/hms/virmerge/files/clusters/vOTUs.tsv")
fpath_vir_merged = here("data/hms/virmerge/files/merge3.tsv")


# ---------- metadata -------------
metadata <- read_tsv(fpath_metadata)
pseq_meta <- metadata %>% 
    column_to_rownames("Project_ID") %>% 
    sample_data()

# ------------ 16S -----------------
amp16s_counts <- read_excel(fpath_16s, sheet = "Counts") %>% 
    column_to_rownames("OTU") %>% 
    as.matrix() %>% 
    otu_table(taxa_are_rows = TRUE)

amp16s_taxa <- read_excel(fpath_16s, sheet = "Taxonomy") %>% 
    column_to_rownames("OTU") %>% 
    as.matrix() %>% 
    tax_table()


pseq_16s <- phyloseq(amp16s_counts, amp16s_taxa, pseq_meta)
save(pseq_16s, file=here("data/processed/pseq_16s.RData"))


# ---------- ITS ---------------

ampITS_counts <- read_excel(fpath_ITS, sheet = "Counts") %>% 
    column_to_rownames("OTU") %>% 
    as.matrix() %>% 
    otu_table(taxa_are_rows = TRUE)

ampITS_taxa <- read_excel(fpath_ITS, sheet = "Taxonomy") %>% 
    column_to_rownames("OTU") %>% 
    as.matrix() %>% 
    tax_table()

pseq_ITS <- phyloseq(ampITS_counts, ampITS_taxa, pseq_meta)
save(pseq_ITS, file=here("data/processed/pseq_ITS.RData"))

# ---------- Virome ------------
# votu <- merge_vOTU_abundance(fpath_vabundance, fpath_vOTUs) %>%
#     column_to_rownames("contig_id_raw") %>%
#     dplyr::select(colnames(.)[colnames(.) %in% metadata$vir_seq_sample_id])

votu <- fread(fpath_vabundance) %>% 
  column_to_rownames("Contig") %>% 
  dplyr::select(colnames(.)[colnames(.) %in% metadata$vir_seq_sample_id])

vmeta <- metadata %>% 
    filter(metadata$vir_seq_sample_id %in% colnames(votu)) %>% 
    column_to_rownames("Project_ID") %>% 
    sample_data()
    
votu <- votu %>% 
    rename(set_names(vmeta$vir_seq_sample_id, rownames(vmeta))) %>% 
    as.matrix() %>% 
    otu_table(taxa_are_rows = TRUE)

vtax <- extract_viral_taxa(fpath_vir_merged) %>%
    column_to_rownames("contig_id_raw") %>% 
    as.matrix() %>% 
    tax_table()

# Create phyloseq object
# pseq_vir <- phyloseq(votu, vtax, pseq_meta)
pseq_vir <- phyloseq(votu, pseq_meta)
save(pseq_vir, file=here("data/processed/pseq_vir.RData"))

