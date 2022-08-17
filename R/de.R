# using "group" variable
edger_de <- function(set, gene_info, ctrs = NULL, normfac = NULL, id = "gene_id", fdr_limit = 0.05, logfc_limit = 0) {
  genes <- gene_info |> 
    select(gene_id, gene_symbol, gene_biotype, description) |> 
    distinct()
  
  cnt <- set$tab[set$sel, ]
  meta <- set$metadata
  groups <- levels(meta$group)
  design_mat <- model.matrix(~0 + group, data = meta)
  colnames(design_mat) <- groups
  
  if (is.null(ctrs)) {
    # all pairwise comparisons
    ctrs <- expand_grid(x = as_factor(groups), y = as_factor(groups)) |>
      filter(as.integer(x) > as.integer(y)) |>
      unite(contrast, c(x, y), sep = "-") |>
      pull(contrast) 
  } 
  contrast_mat <- makeContrasts(contrasts = ctrs, levels = design_mat)
  
  dg <- DGEList(counts = as.matrix(cnt), group = meta$group)
  if (is.null(normfac)) {
    dg <- dg |> calcNormFactors()
  } else {
    n <- tibble(sample = rownames(dg$samples)) |> 
      left_join(normfac, by = c("sample"))
    dg$samples$norm.factors <- n$normfac
  }
    
  fit <- dg |> 
    estimateDisp(design = design_mat) |>
    glmQLFit(design = design_mat)
  
  map_dfr(ctrs, function(ctr) {
    glmQLFTest(fit, contrast = contrast_mat[, ctr]) |>
      topTags(n = 1e16, adjust.method = "BH", sort.by = "none") |>
      pluck("table") |>
      as_tibble(rownames = id) |>
      mutate(contrast = ctr)
  }) |> 
    mutate(contrast = contrast |> str_replace(" - ", "-") |> str_remove_all("group") |> as_factor()) |> 
    left_join(genes, by = "gene_id") |> 
    mutate(
      gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol),
      sig = FDR < fdr_limit & abs(logFC) > logfc_limit
    )
}


# This version of edger_de is with intercept, so it contrasts are against
# baseline. Coefficients and contrasts are automatically calculated, so no
# contrast function necessary.

edger_de_f <- function(set, gene_info, formula, id = "gene_id", fdr_limit = 0.05, logfc_limit = 0) {
  genes <- gene_info |> 
    select(gene_id, gene_symbol, gene_biotype, description) |> 
    distinct()
  
  cnt <- set$tab[set$sel, ]
  meta <- set$metadata
  design_mat <- model.matrix(as.formula(formula), data = meta)
  coefs <- colnames(design_mat)[-1]
  
  fit <- DGEList(counts = as.matrix(cnt)) |>
    calcNormFactors() |>
    estimateDisp(design = design_mat) |>
    glmQLFit(design = design_mat)
  
  map_df(coefs, function(cf) {
    glmQLFTest(fit, coef = cf) |>
      topTags(n = 1e16, adjust.method = "BH", sort.by = "none") |>
      pluck("table") |>
      as_tibble(rownames = id) |>
      mutate(contrast = cf)
  }) |>
    mutate(contrast = factor(contrast, levels = coefs)) |> 
    left_join(genes, by = "gene_id") |> 
    mutate(
      sig = FDR < fdr_limit & abs(logFC) > logfc_limit,
      gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol)
    )
}


de_list <- function(res, group_var, fdr, logfc, fdr_limit = 0.01, logfc_limit = 0, name = "") {
  d <- filter(res, !!sym(fdr) < fdr_limit & abs(!!sym(logfc)) >= logfc_limit) |>
    select(!!group_var, gene_symbol) |>
    group_by(!!sym(group_var))
  kname <- ifelse(name == "", "", paste0(name, ":"))
  ks <- paste0(kname, group_keys(d)[[1]])
  d |>
    distinct() |>
    group_map(~pull(.x, gene_symbol)) |>
    set_names(ks)
}

save_de <- function(de, file) {
  de |> 
    pivot_wider(
      id_cols = c(gene_id, gene_symbol, description),
      names_from = contrast,
      values_from = c(logFC, logCPM, PValue, FDR)
    ) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_tsv(file)
}



# version doing all pairwise contrasts
deseq2_de <- function(set, gene_info, ctrs = NULL, id = "gene_id", fdr_limit = 0.05, logfc_limit = 0) {
  genes <- gene_info |> 
    select(gene_id, gene_symbol, gene_biotype, description) |> 
    distinct()
  
  cnt <- set$tab[set$sel, ]
  meta <- set$metadata
  groups <- levels(meta$group)
  design_mat <- model.matrix(~0 + group, data = meta)
  colnames(design_mat) <- groups
  
  if (is.null(ctrs)) {
    # all pairwise comparisons
    ctrs <- expand_grid(x = as_factor(groups), y = as_factor(groups)) |>
      filter(as.integer(x) > as.integer(y)) |>
      unite(contrast, c(x, y), sep = "-") |>
      pull(contrast) 
  } 

  fit <- DESeqDataSetFromMatrix(countData = as.matrix(cnt), colData = meta, design = ~group) |> 
    DESeq()

  
  map_dfr(ctrs, function(ctr) {
    contr <- str_split(ctr, "-") |> unlist()
    fit |> 
      results(contrast = c("group", contr[1], contr[2])) |> 
      as.data.frame() |> 
      as_tibble(rownames = id) |>
      mutate(contrast = ctr)
  }) |> 
    mutate(contrast = as_factor(contrast)) |> 
    left_join(genes, by = "gene_id") |> 
    rename(logFC = log2FoldChange, PValue = pvalue, FDR = padj) |> 
    mutate(logMean = log2(baseMean)) |> 
    mutate(sig = FDR < fdr_limit & abs(logFC) > logfc_limit)
}




edger_deseq_merge <- function(ed, de) {
  eds <- ed |> 
    select(gene_id, logFC, PValue, FDR, contrast, gene_symbol) |> 
    mutate(tool = "edgeR")
  des <- de |> 
    select(gene_id, logFC, PValue, FDR, contrast, gene_symbol) |> 
    mutate(tool = "DESeq2")
  bind_rows(eds, des)
}

get_rank_k <- function(res, ctrst, direction = 1, max.k) {
  res |>
    filter(contrast == ctrst) |>
    select(gene_id, logFC) |> 
    arrange(desc(direction * logFC)) |> 
    head(max.k) |> 
    mutate(rank  = row_number())
}

plot_compare_top_genes <- function(res1, res2, max.k = 50) {
  map_dfr(levels(res1$contrast), function(ctrst) {
    r1_up <- get_rank_k(res1, ctrst, 1, max.k)
    r2_up <- get_rank_k(res2, ctrst, 1, max.k)
    r1_down <- get_rank_k(res1, ctrst, -1, max.k)
    r2_down <- get_rank_k(res2, ctrst, -1, max.k)
    
    map_dfr(1:max.k, function(k) {
      s1_up <- r1_up |> filter(rank <= k) |> pull(gene_id)
      s2_up <- r2_up |> filter(rank <= k) |> pull(gene_id)
      s1_down <- r1_down |> filter(rank <= k) |> pull(gene_id)
      s2_down <- r2_down |> filter(rank <= k) |> pull(gene_id)
      tibble(
        k = k,
        up = length(intersect(s1_up, s2_up)),
        down = length(intersect(s1_down, s2_down))
      )
    }) |>
      pivot_longer(-k, names_to = "Direction") |> 
      add_column(contrast = ctrst)
  }) |> 
    mutate(contrast = factor(contrast, levels = levels(res1$contrast))) |> 

  ggplot(aes(x = k, y = value, colour = Direction)) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    geom_point() +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max.k + 1)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max.k + 1)) +
    scale_colour_manual(values = okabe_ito_palette) +
    facet_wrap(~contrast, nrow = 1) +
    labs(x = "Gene rank", y = "Genes in common")
}
