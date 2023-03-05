#' Get legend of a plot
#'
#' @param myggplot 
#'
#' @return
#' @export
#'
#' @examples
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Plot richness using phyloseq object
#'
#' @param pseq 
#' @param fout 
#' @param dis_stat 
#'
#' @return
#' @export
#'
#' @examples
plt_richness <- function(pseq, fout, dis_stat=dis_stat) {
  dis_stat <- c("Pre-Tx", "d0", "d7", "d14", "d28", "aGvHD", "Follow Up")
  p1 <- plot_richness(pseq, 
                      measures = c("Chao1"), 
                      x="disease_status",
                      color = "patient_id") +
    scale_x_discrete(limits=dis_stat) +
    stat_summary(fun=abs, colour="black", geom="line", aes(group = 1)) +
    facet_wrap(~patient_id, ncol = 1) +
    theme(legend.position="none", 
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ylab("Chao1") +
    xlab("")
  
  p2 <- plot_richness(pseq, 
                      measures = c("Shannon"), 
                      x="disease_status",
                      color = "patient_id") +
    scale_x_discrete(limits=dis_stat) +
    stat_summary(fun=abs, colour="black", geom="line", aes(group = 1)) +
    facet_wrap(~patient_id, ncol = 1, ) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ylab("Shannon") +
    xlab("")
  
  # legend <- get_legend(p1)
  pdf(fout, width = 5, height = 7)
  p <- grid.arrange(p1, p2, ncol=2, widths=c(2,3), 
                    top = "Chao1 and Shannon diversity", 
                    bottom = "Disease Status")
  print(p)
  dev.off()
}




#' Plot ordination using phyloseq object
#'
#' @param pseq 
#' @param fout 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
plt_ordination <- function(pseq, fout, method) {
  wunifrac_dist = distance(pseq, method=method, weighted=F)
  ordination = ordinate(pseq, method="PCoA", distance=wunifrac_dist)
  
  pdf(fout, width = 6, height = 6)
  p = plot_ordination(pseq, ordination, color="disease_status") + 
    theme(aspect.ratio=1)
  print(p)
  dev.off()
}

plt_taxa_abundance <- function(pseq, taxa_level, fout, dis_stat=dis_stat) {
  pseq_taxa = tax_glom(pseq, taxrank=taxa_level, NArm=FALSE)
  
  pdf(fout, width = 8, height = 6)
  p <- plot_bar(pseq_taxa, fill=taxa_level) +
    facet_grid(~factor(disease_status, levels=dis_stat),
               scale="free_x",
               space = "free_x") +
    theme(strip.text.x = element_text(size = 4))
  print(p)
  dev.off()
}


#' Plot MOFA cum variance
#'
#' @param MOFAobject 
#' @param fout 
#'
#' @return
#' @export
#'
#' @examples
plt_cum_var <- function(MOFAobject, fout) {
  r2 <- MOFAobject@cache$variance_explained$r2_per_factor[[1]]
  r2.dt <- r2 %>%
    as.data.table %>% 
    .[,factor:=as.factor(1:MOFAobject@dimensions$K)] %>%
    melt(id.vars=c("factor"), variable.name="view", value.name = "r2") %>%
    .[,cum_r2:=cumsum(r2), by="view"]
  
  pdf(fout, width = 6, height = 3.5)
  p <- ggline(r2.dt, x="factor", y="cum_r2", color="view") +
    labs(x="Factor number", y="Cumulative variance explained (%)") + 
    theme(legend.title = element_blank(), 
          legend.position = "top",
          axis.text = element_text(size=rel(1)))
  print(p)
  dev.off()
}

#' Plot MOFA weights
#'
#' @param mofa 
#' @param factor 
#' @param view 
#' @param nfeatures 
#'
#' @return
#' @export
#'
#' @examples
plot_weights_fn <- function(mofa, factor=1, view=1, nfeatures=10) {
  p1 <- plot_weights(mofa, 
                     factors = factor, 
                     view = view,
                     nfeatures = nfeatures,
                     text_size = 4
  )
  
  p2 <- plot_top_weights(mofa, 
                         factors = factor, 
                         view = view,
                         nfeatures = nfeatures
  )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  return(p)
}


#' Merge vOTU abundance 
#'
#' @param abundance_contigs_count 
#' @param vOTUs 
#'
#' @return
#' @export
#'
#' @examples
merge_vOTU_abundance <- function(abundance_contigs_count, vOTUs) {
  abundance_all_contigs <- read_table(abundance_contigs_count) %>% 
    dplyr::rename(ctgid=Contig)
  
  votu_abundance <- read_table(vOTUs) %>% 
    left_join(abundance_all_contigs, by = 'ctgid') %>% 
    select(-ctgid) %>% 
    group_by(repid) %>% 
    summarise(across(everything(), sum)) %>% 
    dplyr::rename(contig_id_raw=repid) %>%
    
    return(votu_abundance)
}


#' Extract viral taxonomy
#'
#' @param tbl_merged 
#'
#' @return
#' @export
#'
#' @examples
extract_viral_taxa <- function(tbl_merged) {
  vir_taxa <- read.csv(tbl_merged, sep = "\t") %>% 
    filter(CheckV_miuvig_quality=="High-quality") %>% 
    dplyr::select(c("contig_id_raw", "Viral_Superkingdom", "Viral_Phylum", "Viral_Class", "Viral_Order", "Viral_Family", "Viral_Genus", "Viral_Species", "bin_id")) %>%
    rename_with(~str_replace(., "Viral_", "")) %>%
    # remove prophage contigs
    distinct(contig_id_raw, .keep_all=TRUE)
  
  return(vir_taxa)
}



#' Add NA and then do log10 transformation
#'
#' @param df 
#' @param sampleTable 
#'
#' @return
#' @export
#'
#' @examples
add_NAsamples_log <- function(df, sampleTable) {
  df <- as.data.frame(df)
  df[sampleTable$sample[!sampleTable$sample %in% colnames(df)]] <- NA
  df <- df %>%
    dplyr::select(sort(names(.))) %>%
    as.matrix()
  df <- log10(df+1)
  return(df)
}


#' Perform CLR transformation
#'
#' @param df 
#' @param sampleTable 
#' @param pseudocount 
#'
#' @return
#' @export
#'
#' @examples
clr_trans <- function(df, sampleTable, pseudocount=0) {
  samples_intersection <- intersect(colnames(df), sampleTable$sample)
  df <- df %>%
    dplyr::select(samples_intersection) %>% 
    as.matrix() %>% 
    t()
  df <- df + pseudocount
  df <- df %>% 
    compositions::clr() %>% 
    t() %>% 
    as.data.frame()
  df[sampleTable$sample[!sampleTable$sample %in% colnames(df)]] <- NA
  df <- df %>% 
    dplyr::select(sort(names(.))) %>% 
    as.matrix()
  return(df)
}


#' Read metabolism data
#'
#' @param df_metabolites 
#' @param metadata 
#'
#' @return
#' @export
#'
#' @examples
read_metabolites <- function(df_metabolites, metadata) {
  samples_intersection <- intersect(colnames(df_metabolites), metadata$sample)
  df <- df_metabolites %>%
    dplyr::select(samples_intersection)
  
  df[metadata$sample[!metadata$sample %in% colnames(df)]] <- NA
  df <- df %>% 
    dplyr::select(sort(names(.))) %>% 
    as.matrix()
  return(df)
}


#' Extract VC by checkv quality
#'
#' @param df 
#' @param quality 
#'
#' @return
#' @export
#'
#' @examples
extract_vc_by_quality <- function(df, quality) {
  vcs <- df %>% 
    filter(checkv_quality==quality) %>% 
    dplyr::select(VC.Subcluster) %>% 
    distinct() %>% 
    pull(VC.Subcluster)
  return(vcs)
}


#' Update phyloseq metadata
#'
#' @param pseq 
#' @param fp_metadata 
#' @param sample_order 
#'
#' @return
#' @export
#'
#' @examples
update_physeq_metadata <- function(pseq, fp_metadata, sample_order) {
  mdata <- read_excel(fp_metadata)
  pseq@sam_data <- pseq@sam_data %>% 
    data.frame() %>% 
    rownames_to_column("rowname") %>% 
    left_join(mdata, by = "samplename") %>% 
    column_to_rownames("rowname") %>% 
    sample_data()
  
  pseq <- pseq %>% 
    phyloseq::subset_samples(Timepoint %in% sample_order)
  return(pseq)
}


#' Read CoverM abundance file
#'
#' @param fpath CoverM output file name
#' @param fbin2contig Contig binning to contig mapping
#' @param cov Coverage fraction file or not
#'
#' @return A dataframe
#' @export
#'
#' @examples
read_coverm <- function(fpath, cov = 0, fbin2contig = 0) {
  . <- NULL
  df_abundance <- fread(fpath) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "ds10Ms", "Sample_")) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "Sample_0", "Sample_")) %>%
    setnames("Contig", "genome_id")
  
  # check if bin2contig file is provided
  if (fbin2contig!=0) {
    df_bin2contig <- fread(fbin2contig, header = FALSE, col.names = c("bin_id", "genome_id"))
    df_abundance <- df_abundance %>%
      dplyr::left_join(df_bin2contig, by = "genome_id") %>%
      dplyr::mutate(bin_id = ifelse(is.na(.data$bin_id), .data$genome_id, .data$bin_id)) %>%
      dplyr::mutate(genome_id = .data$bin_id) %>%
      dplyr::select(-.data$bin_id) %>%
      dplyr::group_by(.data$genome_id)
    
    # If read coverage fraction, then choose the maximum value for each bin
    if (cov!=0) {
      df_abundance <- df_abundance %>% summarise_all(max)
    } else {
      df_abundance <- df_abundance %>% summarise_all(sum)
    }
  }
  return(df_abundance)
}


#' Parse vConTACT2 results
#'
#' @param fpath Path of vConTACT2 genome cluster results
#' @param assembler_label "_NODE_" for SPAdes
#'
#' @return list of results
#' @export
#'
#' @examples
read_vConTACT2 <- function(fpath, assembler_label="_NODE_") {
  # TODO: need to tested using a test dataset
  df_vcontact <- fread(fpath) %>%
    mutate(source=ifelse(str_detect(Genome, assembler_label), "queryseq", "refseq")) %>%
    mutate(cluster_status=ifelse(str_detect(`VC Status`, "Overlap"), "Overlap", `VC Status`)) %>%
    mutate(cluster_status=factor(cluster_status)) %>%
    setnames(colnames(.), paste0("vConTACT_", str_replace_all(colnames(.), " ", "_")))
  
  vc2_refs <- df_vcontact %>%
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label, negate = TRUE)) %>%
    dplyr::filter(vConTACT_VC != "")
  
  # Get VC stats
  vcontact_stats <- df_vcontact %>%
    group_by(vConTACT_source, vConTACT_cluster_status) %>%
    dplyr::summarise(seqs_with_VC = sum(vConTACT_VC != ""),
                     seqs_without_VC = sum(vConTACT_VC == "")) %>%
    tidyr::gather(seq_status, num_seqs, -vConTACT_source, -vConTACT_cluster_status)
  
  plot_vc_stats <- ggplot(vcontact_stats, aes(y=num_seqs, x=seq_status, fill=vConTACT_source, label=num_seqs)) +
    geom_point(aes(color=vConTACT_source), alpha=0.5, size=3) +
    facet_grid(vConTACT_cluster_status ~ .) +
    # theme(text = element_text(size = 12)) +
    geom_text_repel()
  
  # Annotate clusters using reference genomes
  vc2_contigs_vclst_anno <- df_vcontact %>%
    # Get distinct cluster ID of contigs in samples
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    dplyr::filter(vConTACT_VC != "") %>%
    dplyr::select(vConTACT_VC) %>%
    distinct() %>%
    # Annotate cluster with reference taxonomy
    inner_join(vc2_refs, by = "vConTACT_VC")
  
  # Whether contig clusters were annotated by reference (1) or not (0)
  df_vcontact <- df_vcontact %>%
    mutate(vConTACT_classified=ifelse(vConTACT_VC %in% vc2_contigs_vclst_anno$vConTACT_VC, 1, 0)) %>%
    # only choose assemblies
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    mutate(vConTACT_VC2=ifelse(vConTACT_VC_Status %in% c("Outlier", "Singleton"), vConTACT_Genome, vConTACT_VC)) %>%
    mutate(vConTACT_VC2=ifelse(str_detect(vConTACT_VC_Status, "Overlap"), vConTACT_VC_Status, vConTACT_VC2)) %>%
    setnames("vConTACT_Genome", "virsorter2_contig_id")
  
  vc2 <- list("vc_tbl" = df_vcontact,
              "vc_stats" = vcontact_stats,
              "vc_plot" = plot_vc_stats,
              "vc_annotated"= vc2_contigs_vclst_anno,
              "vc_refs" = vc2_refs)
  return(vc2)
}

