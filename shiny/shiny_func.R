sh_plot_volcano <- function(d, alpha = 0.05, title = NULL) {
  sres <- d |> filter(FDR <= alpha)
  g <- ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(size = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P), title = title) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
  if (nrow(sres) > 0) {
    g <- g + geom_hline(yintercept = -log10(max(sres$PValue)), colour = "red", linetype = "dashed", alpha = 0.2)
  }
  g
}

sh_plot_ma <- function(d, alpha = 0.05, title = NULL) {
  ggplot(d, aes(x, y)) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = 18)
    ) +
    geom_point(data = d[d$FDR > 0.05,], size = 0.1, colour = "grey50") +
    geom_point(data = d[d$FDR <= 0.05,], size = 0.2, colour = "black") +
    geom_hline(yintercept = 0, colour = "grey70") +
    labs(x = expression(log[2]~CPM), y = expression(log[2]~FC), title = title) 
}




sh_plot_tfe <- function(d) {
  brkn <- c(0, 0.5, 0.8, 0.9, 0.95, 0.99)
  brkn <- c(-brkn, brkn[2:length(brkn)]) |> sort()
  brks <- atanh(brkn)
  ggplot(d, aes(x = x, y = y)) +
    theme_classic() +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_vline(xintercept = 0, colour = "grey50") +
    scale_x_continuous(breaks = brks, labels = brkn) +
    scale_y_continuous(breaks = brks, labels = brkn) +
    labs(x = "Correlation with Tfe1", y = "Correlation with Tfe2")
}


sh_plot_genes <- function(set, genes, y_scale, max_points = 100, what = "count_norm") {
  dat <- set$dat |> 
    filter(gene_id %in% genes)
  n <- nrow(dat)
  if (n == 0) return(NULL)
  
  dat <- dat |> 
    mutate(val = get(what)) |> 
    left_join(set$metadata, by = "sample")
  
  if (y_scale == "log") {
    dat <- dat |> 
      mutate(val = log10(val))
  }
  d <- dat |>
    group_by(sample) |> 
    summarise(
      m = mean(val, na.rm = TRUE),
      se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
      group = first(group),
      replicate = first(replicate) |> as_factor(),
      time = first(time),
      strain = first(strain)
    ) |> 
    replace_na(list(se = 0)) |> 
    mutate(lo = m - se, up = m + se) |> 
    mutate(shape = if_else(m == 0, 24, 21))
  
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
    geom_point(data = d, aes(x = time, y = m, fill = replicate, shape = shape), colour = "grey40", position = pd, size = 4) +
    geom_vline(data = tibble(x = seq(0.5, ntp + 0.5, 1)), aes(xintercept = x), colour = "grey80", alpha = 0.5) +
    viridis::scale_fill_viridis(discrete = TRUE, option = "cividis") +
    xlab("Time (h)")
  if (y_scale == "log") {
    g <- g + labs(y = "Log count")
  } else {
    g <- g +
      labs(y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, NA))
  }
  if (n > 1) {
    g <- g + geom_errorbar(data = d, aes(x = time, ymin = lo, ymax = up, group = replicate), colour = "grey70", alpha = 0.5, position = pd, width = 0.1)
  }
  g
}







# adds gene info and conditions to main data set
sh_prepare_main_data <- function(data) {
  data$dat |>
    left_join(bm_genes, by = "gene_id") |>
    left_join(meta, by = "sample") |> 
    mutate(group = factor(group, levels = levels(star$metadata$group))) |> 
    select(gene_id, gene_symbol, description, sample, group, replicate, count, count_norm, rpkm) |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
}



sh_read_edger_data <- function(path) {
  bm_genes <- read_rds(file.path(path, "bm_genes.rds"))
  star <- read_rds(file.path(path, "star.rds"))
  de <- read_rds(file.path(path, "de.rds"))
  terms <- read_rds(file.path(path, "fterms.rds"))
  tfe_cor <- read_rds(file.path(path, "tfe_cor.rds")) |> 
    left_join(bm_genes, by = "gene_id")
  
  all_genes <- star$genes$gene_id |> unique()
  
  list(
    all_genes = all_genes,
    gene2name = set_names(star$genes$gene_symbol, star$genes$gene_id),
    star = star,
    de = de,
    terms = terms,
    tfe_cor = tfe_cor
  )
}


sh_read_gsea_data <- function(path) {
  bm_genes <- read_rds(file.path(path, "bm_genes.rds"))
  star <- read_rds(file.path(path, "star.rds"))
  terms <- read_rds(file.path(path, "terms.rds"))
  fg <- read_rds(file.path(path, "fg.rds"))
  de <- read_rds(file.path(path, "de.rds"))
  
  list(
    star = star,
    terms = terms,
    fg = fg,
    de = de
  )
  
}
