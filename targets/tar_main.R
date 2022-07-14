targets_main <- function() {
  
  # biomart annotations
  get_annotations <- list(
    tar_target(mart, useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION)),
    tar_target(gns, biomart_fetch_genes(mart)),
    tar_target(bm_gene_lengths, biomart_fetch_gene_lengths(mart, gns$gene_id)),
    tar_target(bm_genes, gns %>% left_join(select(bm_gene_lengths, -chr), by = "gene_id")),
    tar_target(go_terms, go_fetch_go(GO_SPECIES)),
    # tar_target(re_terms, fetch_reactome(mart, star$genes$gene_id)),  # no reactome genes for yeast
    tar_target(kg_terms, get_kegg(species = KEGG_SPECIES, bm_genes = bm_genes, convert_ncbi = FALSE)),
    tar_target(all_terms, list(go = go_terms, kg = kg_terms))
  )

  # directories and metadata
  setup_experiment <- list(
    tar_target(rnaseq_dirs, make_dirs(PATH)),
    tar_target(metadata, make_metadata("tfe", SAMPLE_FILE, SAMPLE_RENAME))
  )

  # quality control
  qc <- list(
    tar_target(fscreen, parse_fscreens(rnaseq_dirs, metadata, suffix = "_R1_screen.txt")),
    tar_target(qcs, parse_qcs(rnaseq_dirs, metadata, paired = TRUE)),
    tar_target(idxstats, parse_idxstats(rnaseq_dirs, metadata)),
    tar_target(tab_star_log, star$star_log),
    
    tar_target(fig_fscreen, plot_fscreen_map(fscreen)),
    tar_target(fig_fscreen_sample, plot_fscreen_sample(fscreen, EXAMPLE)),
    tar_target(fig_read_qual, plot_qualities(qcs)),
    tar_target(fig_read_qual_clust, plot_cluster_qualities(qcs)),
    # tar_target(fig_chrom_proportion, plot_chrom_proportion(idxstats)),
    
    tar_target(fig_star_log, plot_star_log(star$star_log, metadata)),
    tar_target(fig_star_log_map, plot_star_log_map(star$star_log, metadata)),
    tar_target(fig_map_count, plot_mapped_count(star)),
    tar_target(fig_tfe, plot_tfe(idxstats, metadata, x_var = "strain")),
    tar_target(fig_tfe_raw, plot_tfe(idxstats, metadata, x_var = "strain_id")),
    
    tar_target(example_count_file, one_count_file(rnaseq_dirs, 1)),
    tar_target(fig_star_sense, plot_star_sense(example_count_file))
  )
    
  # read star counts
  read_data <- list(
    tar_target(star, read_and_process_star(rnaseq_dirs, metadata, bm_genes, min.count = 10, fix_names_fun = fix_gene_names))
  )
  
  selections <- list(
    tar_target(gene2name, set_names(star$genes$gene_name, star$genes$gene_id)),
    tar_target(edger_sel, edger %>% filter(contrast %in% CONTRAST_SELECTION))
  )

  # read count properties
  count_properties <- list(
    tar_target(fig_sample_dist, plot_sample_distributions(star, x_breaks = c(0, 2, 4), x_lim = c(-1, 5), ncol = 7, colour_var = "strain")),
    tar_target(fig_kernels, plot_kernels(star)),
    tar_target(png_mean_var, plot_mean_var(star) %>% gs("mean_var", 9, 9)),
    tar_target(fig_distance_mat, plot_distance_matrix(star)),
    tar_target(fig_clustering, plot_clustering(star)),
    tar_target(fig_clustering_raw, plot_clustering(star, sample_var = "raw_sample", colour_var = "raw_group")),
    tar_target(fig_pca, plot_pca(star, colour_var = "strain", shape_var = "time")),
    tar_target(fig_umap, plot_umap(star, n_neighbours = 20, min_dist = 0.1, colour_var = "strain", shape_var = "time"))
  )
  
  # differential expression
  differential_expression <- list(
    # DE pairwise group
    tar_target(edger, edger_de(star, gns, fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT) %>% separate_contrasts(metadata)),
    tar_target(deseq, deseq2_de(star, gns, fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT)),
    tar_target(de_cmp, edger_deseq_merge(edger, deseq)),
    
    # DE full model
    tar_target(edger_f, edger_de_f(star, gns, formula = "~ strain + time", fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT)),
    tar_target(edger_fi, edger_de_f(star, gns, formula = "~ strain * time", fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT)),
    tar_target(de_genes, edger_f %>% filter(sig) %>% pull(gene_id) %>% unique()),
    
    # upset
    tar_target(upset_list_edger_fc1, de_list(edger, "contrast", "FDR", "logFC")),
    tar_target(upset_list_deseq_fc1, de_list(deseq, "contrast", "FDR", "logFC")),
    tar_target(upset_list_edger_vs_deseq, de_list(de_cmp, "tool", "FDR", "logFC")),
    
    # DE figures
    tar_target(png_volcano, plot_volcano_grid(edger, metadata) %>% gs("volcano_grid", 12, 12)),
    tar_target(fig_updown, plot_up_down(edger)),
    tar_target(fig_volcano_time, plot_volcano_time(edger)),
    tar_target(fig_volcano_strain, plot_volcano_strain(edger)),
    
    tar_target(fig_volcano_f, plot_volcano(edger_f)),
    tar_target(fig_ma_f, plot_ma(edger_f)),
    tar_target(fig_pdist_f, plot_pdist(edger_f)),
    tar_target(fig_updown_f, plot_up_down(edger_f)),
    tar_target(fig_de_heat, plot_fc_heatmap(star, id_sel = de_genes, max_fc = 1, with_x_text = TRUE, order_col = FALSE))
  )
  
  set_enrichment <- list(
    tar_target(fg_sel, fgsea_all_terms(edger_sel, all_terms)),
    
    tar_target(fig_fg_example_go_0030476, plot_fgsea_enrichment("GO:0030476", edger_sel %>% filter(contrast == "Tfe2_60-WT_60"), go_terms)),
    tar_target(fig_fg_example_go_0032543, plot_fgsea_enrichment("GO:0032543", edger_sel %>% filter(contrast == "Tfe2_60-WT_60"), go_terms))
  )
  
  prepare_for_shiny <- list(
    tar_target(shiny_edger, shiny_data_edger(star, edger, bm_genes, all_terms)),
    tar_target(fg_list, add_links(fg_sel)),
    tar_target(shiny_gsea, shiny_fgsea(star, edger, bm_genes, fg_list)),
    
    tar_target(save_shiny_edger, write_rds(shiny_edger, "shiny/data_edger.rds", compress = "xz")),
    tar_target(save_shiny_gsea, write_rds(shiny_gsea, "shiny/data_gsea.rds", compress = "xz"))
    
  )
  
  info <- list(
    tar_target(edger_version, packageVersion("edgeR"))
  )
  
  c(
    get_annotations,
    setup_experiment,
    selections,
    qc,
    read_data,
    count_properties,
    differential_expression,
    set_enrichment,
    prepare_for_shiny,
    info
  )
}
