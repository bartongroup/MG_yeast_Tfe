shiny_data_edger <- function(star, edger, bm_genes, all_terms) {
  meta <- star$metadata %>% mutate(replicate = as_factor(replicate))
  list(
    metadata = meta,
    dat = star$dat %>%
      left_join(bm_genes, by = "gene_id") %>%
      left_join(meta, by = "sample") %>% 
      mutate(group = factor(group, levels = levels(star$metadata$group))) %>% 
      select(gene_id, gene_name, description, sample, group, replicate, count, rpkm) %>% 
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)),
    genes = star$genes,
    res = edger %>%
      filter(contrast %in% CONTRAST_SELECTION) %>%
      droplevels() %>% 
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)),
    terms = all_terms
  )
}

shiny_fgsea <- function(star, edger, bm_genes, fg_list) {
  meta <- star$metadata %>% mutate(replicate = as_factor(replicate))
  list(
    metadata = meta,
    dat = star$dat %>%
      left_join(bm_genes, by = "gene_id") %>%
      left_join(meta, by = "sample") %>% 
      select(gene_id, gene_name, description, sample, group, replicate, count, rpkm),
    res = edger %>%
      filter(contrast %in% CONTRAST_SELECTION) %>%
      droplevels() %>% 
      mutate(gene_name = if_else(is.na(gene_name), gene_id, gene_name)),
    fgs = fg_list
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
