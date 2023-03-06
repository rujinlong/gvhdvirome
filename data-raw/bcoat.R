## code to prepare `bcoat` dataset goes here

dpath <- here("~/project/tmp/20220718_BCoAT/")
fanno_name2taxid <- here(dpath, "nf_tmp_name2taxid.tsv")
fanno_taxid2lineage <- here(dpath, "nf_tmp_taxaid2lineage.tsv")
fanno <- here(dpath, "nf_tmp_prot2taxa.tsv")
ftree <- here(dpath, "tree.nwk")
# ffaa <- here(dpath, "MSAtrimmed_onlyID.faa")


df_name2taxid <- read.csv(fanno_name2taxid, sep = "\t", header = F, col.names = c("taxname", "taxid"))
df_taxid2lineage <- read.csv(fanno_taxid2lineage, sep = "\t")
df_anno <- read.csv(fanno, sep = "\t", col.names = c("label", "taxname"), header = F) %>% 
  dplyr::left_join(df_name2taxid, by = "taxname") %>% 
  dplyr::left_join(df_taxid2lineage, by = "taxid") %>% 
  mutate(refs=ifelse(taxname %in% c("BAF3", "BAM6"), 1, 0)) %>% 
  mutate(Species = ifelse(is.na(Species), taxname, Species)) %>% 
  mutate(Superkingdom = ifelse(refs==1, "Viruses", Superkingdom)) %>% 
  replace(is.na(.), "") %>% 
  mutate(Superkingdom=paste0("k__", Superkingdom)) %>% 
  mutate(Phylum=paste0("p__", Phylum)) %>% 
  mutate(Class=paste0("c__", Class)) %>% 
  mutate(Order=paste0("o__", Order)) %>% 
  mutate(Family=paste0("f__", Family)) %>% 
  mutate(Genus=paste0("g__", Genus)) %>%
  # mutate(Species=paste0("s__", Species)) %>% 
  mutate(Taxonomy=paste(Superkingdom, Phylum, Class, Order, Family, sep = ";"))


tree <- read.tree(ftree) %>% 
  as_tibble() %>% 
  dplyr::left_join(df_anno, by = "label") %>% 
  dplyr::mutate(Species = ifelse(.data$Species=="BAF3", "VC-1", .data$Species)) %>% 
  dplyr::mutate(Species = ifelse(.data$Species=="BAM6", "VC-2", .data$Species)) %>% 
  as.treedata()


bcoat <- list("tree" = tree,
              "anno" = df_anno)

usethis::use_data(bcoat, overwrite = TRUE)
