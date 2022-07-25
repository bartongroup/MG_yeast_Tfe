# create a few structures for fast selection in enrichment
prepare_for_enrichment <- function(all_terms, all_genes, universes = c("go", "re", "kg")) {
  map(universes, function(u) {
    term_data <- all_terms[[u]]
    if (!is.null(term_data)) {
      # Check for missing term descriptions
      mis_term <- setdiff(term_data$gene2term$term_id, term_data$terms$term_id)
      if (length(mis_term) > 0) {
        dummy <- tibble(
          term_id = mis_term,
          term_name = rep(NA_character_, length(mis_term))
        )
        term_data$terms <- 
          bind_rows(term_data$terms, dummy)
      }
      
      # List to select term info
      term_info <- term_data$terms %>%
        rowwise() %>%
        group_split() %>%
        set_names(term_data$terms$term_id)
      
      # gene-term tibble
      gene_term <- term_data$gene2term %>% 
        # mutate(gene_id = toupper(gene_id)) %>% 
        filter(gene_id %in% all_genes)
      
      # Another named list for selection
      term2gene <- gene_term %>%
        group_by(term_id) %>%
        summarise(genes = list(gene_id)) %>%
        deframe()
      
      list(
        term_info = term_info,
        gene_term = gene_term,
        term2gene = term2gene
      )
    }
  }) %>% 
    set_names(universes)
}



write_rds_name <- function(obj) {
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  file_name <- file.path(path, str_glue("{obj_name}.rds"))
  write_rds(obj, file_name, compress = "xz")
}

save_data_for_shiny <- function(bm_genes, star, edger_sel, edger_fi, tfe_cor,
                                all_terms, fg_sel, fg_fi, fg_tfe) {
  de <- list(sel = edger_sel, fi = edger_fi)
  fg <- list(sel = fg_sel, fi = fg_fi, tfe = fg_tfe)
  all_genes <- star$genes$gene_id %>% unique()
  terms <- prepare_for_enrichment(all_terms, all_genes)
  tfe_cor <- tfe_cor %>%
    pivot_wider(id_cols = gene_id, names_from = contrast, values_from = correlation)
  
  write_rds_name(bm_genes)
  write_rds_name(star)
  write_rds_name(de)
  write_rds_name(tfe_cor)
  write_rds_name(terms)
  write_rds_name(fg)
}


shiny_data_edger <- function(star, edger_sel, edger_fi, tfe_cor, bm_genes, all_terms) {
  meta <- star$metadata %>% mutate(replicate = as_factor(replicate))
  list(
    metadata = meta,
    gene_info = bm_genes,
    dat = star$dat %>%
      left_join(bm_genes, by = "gene_id") %>%
      left_join(meta, by = "sample") %>% 
      mutate(group = factor(group, levels = levels(star$metadata$group))) %>% 
      select(gene_id, gene_name, description, sample, group, replicate, count, count_norm, rpkm) %>% 
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)),
    genes = star$genes,
    de = edger %>%
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)) %>% 
      select(contrast, gene_id, gene_name, logFC, logCPM, PValue, FDR, gene_biotype, description),
    tfe = tfe_cor %>% 
      pivot_wider(id_cols = gene_id, names_from = "contrast", values_from = "correlation") %>% 
      left_join(bm_genes, by = "gene_id") %>% 
      select(gene_id, gene_name, gene_biotype, description, corTfe1, corTfe2),
    terms = all_terms
  )
}

shiny_fgsea <- function(star, edger, bm_genes, fg_list, fg_tfe) {
  stopifnot(names(fg_list) == names(fg_tfe))
  meta <- star$metadata %>% mutate(replicate = as_factor(replicate))
  univs <- names(fg_list)
  list(
    metadata = meta,
    dat = star$dat %>%
      left_join(bm_genes, by = "gene_id") %>%
      left_join(meta, by = "sample") %>% 
      select(gene_id, gene_name, description, sample, group, replicate, count, rpkm),
    de = edger %>%
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)),
    fgs = map(univs, function(u) {
      bind_rows(fg_list[[u]], fg_tfe[[u]])
    }) %>% set_names(univs)
  )
}






add_links <- function(fg_lst) {
  go_url <- "http://amigo.geneontology.org/amigo/term/{term}"
  re_url <- "https://reactome.org/content/detail/{term}"
  kg_url <- "https://www.genome.jp/pathway/{term}"
  
  list(
    go = fg_lst$go %>% mutate(term = glue::glue('<a href=http://amigo.geneontology.org/amigo/term/{term}>{term}</a>')),
    # gs = fg_lst$gs %>% mutate(term = glue::glue('<a href=http://amigo.geneontology.org/amigo/term/{term}>{term}</a>')),
    re = fg_lst$re %>% mutate(term = glue::glue('<a href=https://reactome.org/content/detail/{term}>{term}</a>')),
    kg = fg_lst$kg %>% mutate(term = glue::glue('<a href=https://www.genome.jp/pathway/{term}>{term}</a>'))
  )
}
