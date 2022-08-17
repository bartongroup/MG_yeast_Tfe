# Spike-in normalisation. For each spike-in, spike_rat is the ratio of the count
# to the mean count in the spike-in. Then, as we have three spike-ins, a
# normalising factor for each sample is calculated as the mean of the three
# spike_rats.

spike_norm <- function(set, spikes, what="count") {
  set$dat |> 
    filter(gene_id %in% spikes) |> 
    mutate(value = get(what)) |> 
    select(gene_id, sample, value) |> 
    group_by(gene_id) |> 
    mutate(spike_rat = value / mean(value)) |> 
    ungroup() |> 
    group_by(sample) |> 
    mutate(normfac = mean(spike_rat)) |> 
    ungroup()
}


normalise_by_spike <- function(set, spikes) {
  sp <- spike_norm(set, spikes)
  nrm <- sp |> 
    select(sample, normfac) |> 
    distinct()
  
  set$spike_normfac <- nrm
  set$dat <- set$dat |> 
    left_join(nrm, by = "sample") |> 
    mutate(count_spikenorm = count / normfac) |> 
    select(-normfac)
  set$spikeins <- spikes
  
  set
  
}

plot_spike_norm <- function(set, spikes=SPIKE_IDS) {
  sp <- spike_norm(set, spikes)
  sp |> 
    ggplot(aes(x=sample, y=spike_rat, colour=gene_id)) +
    geom_point() +
    geom_point(aes(y=normfac), shape=3, size=3) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5) 
    ) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x=NULL, y="Spike ratio")
}


plot_spike_comparison <- function(set, nrow=1) {
  spikes <- set$spikeins
  dat <- set$dat |> 
    filter(gene_id %in% spikes) |> 
    left_join(set$metadata, by = "sample") |> 
    select(gene_id, sample, count, group) |> 
    arrange(sample)
  
  # splitting into all combinations of x and y for a pair plot
  expand_grid(x = spikes, y = spikes) |> 
    filter(x < y) |> 
    unite(name, c(x, y), sep=":", remove=FALSE) |> 
    rowwise() |> 
    mutate(
      val_x = list(dat |> filter(gene_id == x) |> select(count_x = count, group_x = group)),
      val_y = list(dat |> filter(gene_id == y) |> select(count_y = count, group_y = group))
    ) |> 
    unnest(c(val_x, val_y)) |> 
  ggplot(aes(x = count_x, y = count_y, colour = group_x)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_point() +
    scale_colour_manual(values = okabe_ito_palette) +
    facet_wrap(~name, nrow=nrow) +
    labs(x=NULL, y=NULL, colour="Group")
}


plot_spike_cv <- function(set, min.m = 10) {
  d <- set$dat |>
    filter(good) |> 
    left_join(select(set$metadata, sample, group), by="sample") |>
    mutate(gene_id = as_factor(gene_id)) |> 
    group_by(gene_id, group) |>
    summarise(
      M_cn = mean(count_norm), S_cn = sd(count_norm), CV_cn = S_cn / M_cn,
      M_sp = mean(count_spikenorm), S_sp = sd(count_spikenorm), CV_sp = S_sp / M_sp,
    ) |> 
    filter(M_cn > min.m & M_sp > min.m) |> 
    ungroup()
  dm <- d |> 
    group_by(group) |> 
    summarise(perc = 100 * length(which(CV_sp > CV_cn)) / n()) |> 
    mutate(perc = sprintf("%4.1f", perc))
  d |> 
    ggplot(aes(x = CV_cn, y = CV_sp)) +
    theme_bw() +
    geom_point(alpha=0.1) +
    geom_abline(slope=1, intercept=0, colour="red") +
    geom_text(x = 0.3, y=1.85, aes(label = perc), data = dm) +
    facet_wrap(~group) +
    labs(x = "CV normalised to library size", y = "CV normalised by spike-ins")
}


plot_spike_ratios <- function(set) {
  spikes <- set$spikeins
  dat <- set$dat |> 
    filter(gene_id %in% spikes) |> 
    left_join(set$metadata, by = "sample") |> 
    select(gene_id, sample, count, group) |> 
    arrange(sample)
  
  
  expand_grid(x = spikes, y = spikes) |> 
    filter(x < y) |> 
    unite(name, c(y, x), sep=":", remove=FALSE) |> 
    rowwise() |> 
    mutate(
      val_x = list(dat |> filter(gene_id == x) |> select(count_x = count, sample_x = sample)),
      val_y = list(dat |> filter(gene_id == y) |> select(count_y = count, sample_y = sample))
    ) |> 
    unnest(c(val_x, val_y)) |> 
    mutate(r = count_y / count_x, sample = factor(sample_x, levels = set$metadata$sample)) |> 
  ggplot(aes(x = sample, y = r, colour = name, group = name)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5) 
    ) +
    geom_line() +
    geom_point() +
    scale_colour_brewer(type = "qual", palette = "Dark2") +
    labs(x = NULL, y = "Spike-to-spike ratio", colour = "Spike pair")
}


plot_spike_fraction <- function(set, norm = "input") {
  spikes <- set$spikeins
  
  if(norm == "input") {
    libsize <- set$input_normfac
  } else {
    libsize <- set$mapped_normfac
  }

  set$dat |> 
    filter(gene_id %in% spikes) |> 
    left_join(set$metadata, by = "sample") |> 
    select(gene_id, sample, count, group) |> 
    left_join(libsize, by="sample") |> 
    mutate(frac = 100 * count / size) |> 
  ggplot(aes(x = sample, y = frac, colour = gene_id, group = gene_id)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5) 
    ) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values = okabe_ito_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    labs(x = NULL, y = "Percetnage of total reads", colour = "Spike")

}
