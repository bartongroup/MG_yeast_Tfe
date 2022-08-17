separate_contrasts <- function(res, meta) {
  mt <- meta |> 
    select(group, strain, time) |> 
    distinct()
  
  res |> 
    separate(contrast, c("group1", "group2"), remove = FALSE, sep = "-") |> 
    left_join(mt, by = c("group1" = "group")) |> rename(time1 = time, strain1 = strain) |> 
    left_join(mt, by = c("group2" = "group")) |> rename(time2 = time, strain2 = strain) 
}



plot_v_grid <- function(r) {
  ggplot(r, aes(x = x, y = y, colour = sig)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_colour_manual(values = c("grey70", "black")) +
    facet_grid(y_var ~ x_var) +
    geom_vline(xintercept = 0, size = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}



plot_volcano_grid <- function(res, meta) {
  groups <- levels(meta$group)
  r <- res |> 
    mutate(
      group1 = factor(group1, levels = groups),
      group2 = factor(group2, levels = groups)
    ) |> 
    mutate(
      x = logFC,
      y = -log10(PValue)
    ) |> 
    select(x, y, sig, x_var = group1, y_var = group2)
  rm(res)
  plot_v_grid(r)
}

plot_volcano_time <- function(res) {
  r <- res |> 
    filter(time1 == time2) |> 
    mutate(
      contrast = str_remove_all(contrast, "_\\d+"),
      x = logFC,
      y = -log10(PValue)
    ) |> 
    select(x, y, sig, x_var = contrast, y_var = time1)
  rm(res)
  plot_v_grid(r)
}


plot_volcano_strain <- function(res) {
  r <- res |> 
    filter(strain1 == strain2) |> 
    mutate(
      contrast = str_remove_all(contrast, "\\w+_"),
      x = logFC,
      y = -log10(PValue)
    ) |> 
    select(x, y, sig, x_var = contrast, y_var = strain1)
  rm(res)
  plot_v_grid(r)
}


plot_tfe <- function(idxstats, meta, x_var = "strain") {
  idxstats |> 
    filter(chr %in% c("Tfe1", "Tfe2")) |> 
    left_join(meta, by = "raw_sample") |> 
    unite(timerep, c(time, replicate), remove = FALSE) |> 
  ggplot(aes(x = timerep, y = count, colour = time)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_segment(aes(xend = timerep, yend = 0), size = 1.1, alpha = 0.5) +
    geom_point() +
    facet_grid(as.formula(str_glue("chr ~ {x_var}"))) +
    scale_colour_manual(values = okabe_ito_palette) +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 5, 10, 20, 50, 100) * 1000, labels = scales::comma) +
    labs(x = "Time / replicate", y = "Normalised cound per Tfe reference")
}


normalise_to_wt <- function(set) {
  nrm <- set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    filter(strain == "WT") |> 
    group_by(gene_id, time) |> 
    summarise(mean_wt = mean(count_norm) + 0.5) |> 
    ungroup(time) |> 
    mutate(normfac = mean_wt / mean(mean_wt)) |> 
    ungroup()
  set$dat <- set$dat |> 
    left_join(set$metadata |> select(sample, time), by = "sample") |> 
    left_join(nrm, by = c("gene_id", "time")) |> 
    mutate(count_wt = count_norm / normfac, .after = count_norm) |> 
    select(-c(time, mean_wt, normfac))
  set
}


plot_gene_time <- function(set, genes, val = "count_norm", max_points = 100) {
  dat <- set$dat |> 
    filter(gene_id %in% genes) 
  n <- nrow(dat)
  if (n == 0) return(NULL)
  
  d <- dat |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(value = get(val))
  
  pd <- ggplot2::position_dodge(width = 0.4)
  
  th <- theme_bw() + theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "none"
  )
  
  
  ntp <- length(unique(d$time))
  g <- ggplot() +
    th +
    scale_shape_identity() +  # necessary for shape mapping
    facet_grid(strain ~ .) +
    geom_point(data = d, aes(x = time, y = value, fill = replicate, shape = 21), colour = "grey40", position = pd, size = 4) +
    geom_vline(data = tibble(x = seq(0.5, ntp + 0.5, 1)), aes(xintercept = x), colour = "grey80", alpha = 0.5) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "cividis") +
    labs(x = "Time (h)", y = val) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, NA))
  g
}

get_tfe <- function(set, idxs) {
  w <- idxs |> 
    filter(chr %in% c("Tfe1", "Tfe2")) |> 
    left_join(set$mapped_normfac, by = "sample") |> 
    mutate(count_norm = count / normfac) 
  w |> 
    pivot_wider(id_cols = sample, names_from = chr, values_from = count) |> 
    column_to_rownames("sample") |>
    as.matrix() |> 
    DESeq2::rlog() |>
    as_tibble(rownames = "sample") |> 
    pivot_longer(-sample, names_to = "strain", values_to = "rlog") |> 
    left_join(select(w, sample, strain = chr, count, count_norm), by = c("sample", "strain"))
}


tfe_correlation <- function(set, tfe, what = "count_norm") {
  tfe_tab <- tfe |>
    mutate(value = get(what)) |> 
    pivot_wider(id_cols = sample, names_from = strain, values_from = value) |> 
    column_to_rownames("sample") |>
    as.matrix()
  data_tab <- dat2mat(set$dat, what)
  cor(t(data_tab), tfe_tab) |> 
    as_tibble(rownames = "gene_id") |> 
    pivot_longer(-gene_id, names_to = "contrast", values_to = "correlation") |> 
    mutate(contrast = recode(contrast, Tfe1 = "corTfe1", Tfe2 = "corTfe2"))
}


plot_tfe_correlation <- function(tfe_cor, tr = "atanh") {
  brks <- c(0, 0.5, 0.8, 0.9, 0.95, 0.99)
  brks <- c(-brks, brks)
  tfe_cor |> 
    pivot_wider(id_cols = gene_id, names_from = contrast, values_from = correlation) |> 
  ggplot(aes(x = corTfe1, y = corTfe2)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    scale_x_continuous(trans = tr, breaks = brks) +
    scale_y_continuous(trans = tr, breaks = brks) +
    labs(x = "Correlation with Tfe1", y = "Correlation with Tfe2")
}

plot_gene_tfe <- function(set, tfe, gid, what = "count_norm") {
  set$dat |> 
    filter(gene_id == gid) |> 
    left_join(set$genes, by = "gene_id") |> 
    select(sample, strain = gene_symbol, rlog, count, count_norm) |> 
    bind_rows(tfe) |> 
    mutate(value = get(what)) |> 
    mutate(
      sample = factor(sample, levels = set$metadata$sample),
      strain = as_factor(strain)
    ) |> 
  ggplot(aes(x = sample, y = value)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_hline(yintercept = 0, colour = "grey50")  +
    geom_segment(aes(xend = sample, yend = 0), colour = "grey80") +
    geom_point() +
    facet_wrap(~ strain, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = what)
}

make_tf_cor_table <- function(tfe_cor, bm_genes) {
  tfe_cor |> 
    mutate(correlation = signif(correlation, 3)) |> 
    pivot_wider(id_cols = gene_id, names_from = contrast, values_from = correlation) |> 
    left_join(bm_genes |> select(gene_id, gene_symbol, gene_biotype, description))
}



get_sgd_data <- function(gids) {
  api <- "https://www.yeastgenome.org/backend/locus/"
  map_dfr(gids, function(gid) {
    res <- httr::GET(paste0(api, gid))
    d <- jsonlite::fromJSON(rawToChar(res$content))
    c(
      gene_id = gid,
      gene_symbol = d$gene_name,
      description = d$name_description,
      description_full = d$description
    )
  })
}