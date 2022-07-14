separate_contrasts <- function(res, meta) {
  mt <- meta %>% 
    select(group, strain, time) %>% 
    distinct()
  
  res %>% 
    separate(contrast, c("group1", "group2"), remove = FALSE, sep = "-") %>% 
    left_join(mt, by = c("group1" = "group")) %>% rename(time1 = time, strain1 = strain) %>% 
    left_join(mt, by = c("group2" = "group")) %>% rename(time2 = time, strain2 = strain) 
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
  r <- res %>% 
    mutate(
      group1 = factor(group1, levels = groups),
      group2 = factor(group2, levels = groups)
    ) %>% 
    mutate(
      x = logFC,
      y = -log10(PValue)
    ) %>% 
    select(x, y, sig, x_var = group1, y_var = group2)
  rm(res)
  plot_v_grid(r)
}

plot_volcano_time <- function(res) {
  r <- res %>% 
    filter(time1 == time2) %>% 
    mutate(
      contrast = str_remove_all(contrast, "_\\d+"),
      x = logFC,
      y = -log10(PValue)
    ) %>% 
    select(x, y, sig, x_var = contrast, y_var = time1)
  rm(res)
  plot_v_grid(r)
}


plot_volcano_strain <- function(res) {
  r <- res %>% 
    filter(strain1 == strain2) %>% 
    mutate(
      contrast = str_remove_all(contrast, "\\w+_"),
      x = logFC,
      y = -log10(PValue)
    ) %>% 
    select(x, y, sig, x_var = contrast, y_var = strain1)
  rm(res)
  plot_v_grid(r)
}


plot_tfe <- function(idxstats, meta, x_var = "strain") {
  idxstats %>% 
    filter(chr %in% c("Tfe1", "Tfe2")) %>% 
    left_join(meta, by = "raw_sample") %>% 
    unite(timerep, c(time, replicate), remove = FALSE) %>% 
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
  nrm <- set$dat %>% 
    left_join(set$metadata, by = "sample") %>% 
    filter(strain == "WT") %>% 
    group_by(gene_id, time) %>% 
    summarise(mean_wt = mean(count_norm) + 0.5) %>% 
    ungroup(time) %>% 
    mutate(normfac = mean_wt / mean(mean_wt)) %>% 
    ungroup()
  set$dat <- set$dat %>% 
    left_join(set$metadata %>% select(sample, time), by = "sample") %>% 
    left_join(nrm, by = c("gene_id", "time")) %>% 
    mutate(count_wt = count_norm / normfac, .after = count_norm) %>% 
    select(-c(time, mean_wt, normfac))
  set
}


plot_gene_time <- function(set, genes, val = "count_norm", max_points = 100) {
  dat <- set$dat %>% 
    filter(gene_id %in% genes) 
  n <- nrow(dat)
  if (n == 0) return(NULL)
  
  d <- dat %>% 
    left_join(set$metadata, by = "sample") %>% 
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
