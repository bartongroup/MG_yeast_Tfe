write_rds_name <- function(obj) {
  path <- file.path("shiny", "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  file_name <- file.path(path, str_glue("{obj_name}.rds"))
  write_rds(obj, file_name, compress = "xz")
}

save_data_for_shiny <- function(bm_genes, star, edger_sel, edger_fi, tfe_cor,
                                fterms, fg_sel, fg_fi, fg_tfe) {
  de <- list(sel = edger_sel, fi = edger_fi)
  fg <- list(sel = fg_sel, fi = fg_fi, tfe = fg_tfe)
  all_genes <- star$genes$gene_id |> unique()
  terms <- fterms
  tfe_cor <- tfe_cor |>
    pivot_wider(id_cols = gene_id, names_from = contrast, values_from = correlation)
  
  write_rds_name(bm_genes)
  write_rds_name(star)
  write_rds_name(de)
  write_rds_name(tfe_cor)
  write_rds_name(terms)
  write_rds_name(fg)
}


shiny_data_edger <- function(star, edger_sel, edger_fi, tfe_cor, bm_genes, all_terms) {
  meta <- star$metadata |> mutate(replicate = as_factor(replicate))
  list(
    metadata = meta,
    gene_info = bm_genes,
    dat = star$dat |>
      left_join(bm_genes, by = "gene_id") |>
      left_join(meta, by = "sample") |> 
      mutate(group = factor(group, levels = levels(star$metadata$group))) |> 
      select(gene_id, gene_symbol, description, sample, group, replicate, count, count_norm, rpkm) |> 
      mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol)),
    genes = star$genes,
    de = edger |>
      mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol)) |> 
      select(contrast, gene_id, gene_symbol, logFC, logCPM, PValue, FDR, gene_biotype, description),
    tfe = tfe_cor |> 
      pivot_wider(id_cols = gene_id, names_from = "contrast", values_from = "correlation") |> 
      left_join(bm_genes, by = "gene_id") |> 
      select(gene_id, gene_symbol, gene_biotype, description, corTfe1, corTfe2),
    terms = all_terms
  )
}

shiny_fgsea <- function(star, edger, bm_genes, fg_list, fg_tfe) {
  stopifnot(names(fg_list) == names(fg_tfe))
  meta <- star$metadata |> mutate(replicate = as_factor(replicate))
  univs <- names(fg_list)
  list(
    metadata = meta,
    dat = star$dat |>
      left_join(bm_genes, by = "gene_id") |>
      left_join(meta, by = "sample") |> 
      select(gene_id, gene_symbol, description, sample, group, replicate, count, rpkm),
    de = edger |>
      mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol)),
    fgs = map(univs, function(u) {
      bind_rows(fg_list[[u]], fg_tfe[[u]])
    }) |> set_names(univs)
  )
}






add_links <- function(fg_lst) {
  list(
    go = fg_lst$go |> mutate(term_id = str_glue('<a href=http://amigo.geneontology.org/amigo/term/{term_id}>{term_id}</a>')),
    re = fg_lst$re |> mutate(term_id = str_glue('<a href=https://reactome.org/content/detail/{term_id}>{term_id}</a>')),
    kg = fg_lst$kg |> mutate(term_id = str_glue('<a href=https://www.genome.jp/pathway/{term_id}>{term_id}</a>')),
    bp = fg_lst$bp |> mutate(term_id = str_glue("<a href=https://tripod.nih.gov/bioplanet/detail.jsp?pid=bioplanet_{term_id}&target=pathway>{term_id}</a>"))
  )
}
