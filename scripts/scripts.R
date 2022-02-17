library(tidyverse)
library(readxl)
library(gridExtra)
library(reshape2)
library(data.table)
library(ggpubr)
library(speedyseq)

dis_stat_vir <- c("Pre-Tx", "d0", "d7", "d14", "d28", "aGvHD")
dis_stat_16s <- c("Pre-Tx", "d0", "d7", "d14", "d28", "aGvHD", "Follow Up")
dis_stat_18s <- c("Pre-TX", "d0", "d7", "d14", "d28", "aGvHD", "Follow_up")

get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

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


extract_viral_taxa <- function(tbl_merged) {
    vir_taxa <- read.csv(tbl_merged, sep = "\t") %>% 
        filter(CheckV_miuvig_quality=="High-quality") %>% 
        dplyr::select(c("contig_id_raw", "Viral_Superkingdom", "Viral_Phylum", "Viral_Class", "Viral_Order", "Viral_Family", "Viral_Genus", "Viral_Species", "bin_id")) %>%
        rename_with(~str_replace(., "Viral_", "")) %>%
        # remove prophage contigs
        distinct(contig_id_raw, .keep_all=TRUE)

    return(vir_taxa)
}


# MOFA2 analyses
add_NAsamples_log <- function(df, sampleTable) {
  df <- as.data.frame(df)
  df[sampleTable$sample[!sampleTable$sample %in% colnames(df)]] <- NA
  df <- df %>%
    dplyr::select(sort(names(.))) %>%
    as.matrix()
  df <- log10(df+1)
  return(df)
}

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


extract_vc_by_quality <- function(df, quality) {
  vcs <- df %>% 
    filter(checkv_quality==quality) %>% 
    dplyr::select(VC.Subcluster) %>% 
    distinct() %>% 
    pull(VC.Subcluster)
  return(vcs)
}



# ------ Amplicon --------
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
