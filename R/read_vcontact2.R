#' Title
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
