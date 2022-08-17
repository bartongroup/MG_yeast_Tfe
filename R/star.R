one_count_file <- function(paths, n) {
  path <- paths$readcount
  files <- list.files(path, pattern = ".*txt", full.names = TRUE)
  files[n]
}


read_and_process_star <- function(paths, meta, gene_info, suffix = ".txt", min.count = 10, fix_names_fun = NULL) {
  star.column <- star_strand_column(paths, meta, suffix)
  cat(paste0("Detected strand - ", names(star.column), ". Reading STAR column ", star.column, ".\n"))
  
  set <- parse_star_counts(paths, meta, star.column, suffix, fix_names_fun) |> 
    add_star_log(paths, meta) |> 
    add_gene_names(gene_info) |> 
    normalize_star_counts(gene_info) |> 
    filter_star_min_count(min.count)

  return(set)
}

add_star_log <- function(set, paths, meta) {
  set$star_log <- parse_star_logs(paths, meta)
  return(set)
}

parse_one_star_log <- function(file, smpl) {
  if (!file.exists(file)) {
    warning(paste("File", file, "not found"))
    return(NULL)
  }
  tibble(x = readLines(file)) |> 
    separate(x, c("key", "value"), sep = "\\t", fill = "right") |> 
    drop_na() |> 
    mutate(key = str_remove(key, " \\|")) |> 
    mutate(key = str_remove(key, "^\\s+")) |> 
    mutate(raw_sample = smpl)
}

# parse all star logs
parse_star_logs <- function(paths, meta) {
  path <- paths$starmap
  s2n <- set_names(meta$sample, meta$raw_sample)
  meta$raw_sample |> 
    map_dfr(~parse_one_star_log(file.path(path, paste0(.x, "_Log.final.out")), .x)) |> 
    mutate(sample = as.character(s2n[raw_sample]))
}

# tabulate star log for printing
tabulate_star_log <- function(slog) {
  slog |> 
    select(-raw_sample) |> 
    pivot_wider(names_from = "sample", values_from = "value") |> 
    rename(Description = key) 
}

parse_one_star_count <- function(file, smpl, column = 2, fix_names_fun = NULL) {
  if (!file.exists(file)) {
    warning(str_glue("File {file} not found."))
    return(NULL)
  }
  d <- read_tsv(file, col_names = FALSE, skip = 4, col_types = "ciii") |> 
    select(1, all_of(column)) |> 
    set_names("gene_id", "count") |> 
    add_column(raw_sample = smpl, .after = "gene_id")
  if (!is.null(fix_names_fun)) {
    d <- d |> mutate(gene_id = fix_names_fun(gene_id))
  }
}

# Parse all star count files. Returns a list of tab (wide matrix) and dat (long tibble).
parse_star_counts <- function(paths, meta, column = 2, suffix = ".txt", fix_names_fun) {
  path <- paths$readcount
  s2n <- set_names(meta$sample, meta$raw_sample)
  dat <- meta$raw_sample |> 
    map_dfr(~parse_one_star_count(file.path(path, paste0(.x, suffix)), .x, column, fix_names_fun)) |> 
    mutate(sample = as.character(s2n[raw_sample])) |> 
    group_by(gene_id) |>
    mutate(gene_count = sum(count)) |>
    ungroup() |>
    filter(gene_count > 0) |>
    select(gene_id, sample, count) |> 
    mutate(sample = factor(sample, levels = meta$sample))
  tab <- dat |> 
    pivot_wider(id_cols = "gene_id", names_from = "sample", values_from = "count") |> 
    column_to_rownames("gene_id") |> 
    as.matrix()
  list(dat = dat, tab = tab, metadata = meta, sel = rownames(tab))
}


parse_idxstats <- function(paths, meta) {
  path <- paths$chrcount
  map2_dfr(meta$raw_sample, meta$sample, function(rsam, sam) {
    sfile <- file.path(path, glue::glue("{rsam}.txt"))
    if (file.exists(sfile)) {
    read_tsv(sfile, col_names = c("chr", "length", "count", "unmapped"), show_col_types = FALSE) |> 
      add_column(raw_sample = rsam, sample = sam)
    }
  })
}


add_gene_names <- function(set, gene_info) {
  set$genes <- tibble(gene_id = rownames(set$tab)) |> 
    left_join(gene_info |> select(gene_id, gene_symbol) |> distinct(), by = "gene_id") |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
  set
}

# used in normalise_to_library
normalise_to_size <- function(set, libsize = NULL) {
  if (is.null(libsize)) {
    libsize <- set$dat |> 
      group_by(sample) |> 
      summarise(size = sum(count, na.rm = TRUE)) 
  }
  
  libsize <- libsize |> mutate(normfac = size / mean(size))
  
  # look-up tables are much faster than left_join
  list(
    norm_fac = set_names(libsize$normfac, libsize$sample),
    size_fac = set_names(libsize$size / 1e6, libsize$sample),
    normfac_tab = libsize
  )
  
}

# normalize to library size and RPKM
normalise_to_library <- function(set, gene_info, input_size) {
  
  len_fac = set_names(gene_info$length / 1e3, gene_info$gene_id)
  
  # normalised to total mapped and counted reads
  mapped <- normalise_to_size(set)
  # normalised to total input reads
  input <- normalise_to_size(set, input_size)
  
  
  set$dat <- set$dat |> 
    mutate(
      count_norm = count / mapped$norm_fac[sample] |> unname(),
      rpkm = (count + 1) / (mapped$size_fac[sample] * len_fac[gene_id] |> unname())
      #count_inputnorm = count / input$norm_fac[sample] |> unname,
      #rpkm_inputnorm = (count + 1) / (input$size_fac[sample] * len_fac[gene_id] |> unname),
    )
  
  set$mapped_normfac <- mapped$normfac_tab
  set$input_normfac <- input$normfac_tab
  
  return(set)
}


normalise_edger <- function(set) {
  ed <- edgeR::DGEList(set$tab) |> 
    edgeR::calcNormFactors() |> 
    pluck("samples") |> 
    rownames_to_column("sample") |> 
    rename(normfac = norm.factors) |> 
    select(sample, normfac) |> 
    as_tibble()
  
  # look-up tables are much faster than left_join
  norm_fac <- set_names(ed$normfac, ed$sample)
  
  set$dat <- set$dat |> 
    mutate(count_tmm = count / norm_fac[sample] |> unname())
  set$edger_normfac <- ed
  
  return(set)
}

# regularised log from DESeq2 (originally log2, we convert to log10)
regularised_log <- function(set) {
  rdat <- DESeq2::rlog(set$tab) |> 
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample", values_to = "rlog") |> 
    mutate(
      rlog = rlog / log2(10),
      sample = factor(sample, levels = set$metadata$sample)
    )
  set$dat <- set$dat |> 
    left_join(rdat, by = c("gene_id", "sample"))
 
  return(set) 
}


# edit: added DESeq2's regularised logarithm
normalize_star_counts <- function(set, gene_info) {
  libsize <- set$star_log |> 
    filter(key == "Number of input reads") |> 
    mutate(size = as.numeric(value)) |> 
    select(sample, size)
  
  set |> 
    normalise_to_library(gene_info, libsize) |> 
    #normalise_edger() |> 
    regularised_log() |> 
    normalise_to_wt()
}



plot_star_log <- function(slog, meta, 
                        descs = c(
                          "Number of input reads",
                          "Uniquely mapped reads %",
                          "% of reads mapped to multiple loci",
                          "% of reads unmapped: too short")
                        ) {
  s <- slog |> 
    filter(key %in% descs) |>
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(key = factor(key, levels = descs)) |> 
    select(-raw_sample)
  su <- s |> 
    filter(key == "Uniquely mapped reads %") |> 
    arrange(value)
  
  s |> 
    mutate(sample = factor(sample, levels = su$sample)) |> 
  ggplot(aes(x = sample, y = value, group = key)) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank()) +
    geom_segment(aes(xend = sample, yend = 0), colour = "grey70") +
    geom_point() +
    facet_grid(~key, scales = "free") +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    geom_text(aes(y = value * 1.05, label = "")) # blank geom to expand axis
}

plot_star_log_map <- function(slog, meta) {
  slog |> 
    filter(str_detect(key, "%")) |> 
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(value = na_if(value, 0)) |> 
    mutate(sample = factor(sample, levels = meta$sample)) |> 
    mutate(Description = factor(key)) |> 
  ggplot(aes(x = sample, y = Description, fill = value)) + 
    theme_bw() +
    geom_tile() +
    scale_fill_viridis_c(option = "cividis", trans = "log10", limits = c(0.1, 100), breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(x = NULL, y = NULL, fill = "Percentage")
}


plot_star_sense <- function(file) {
  read_tsv(file, skip = 4, col_names = c("gene", "Unstranded", "First", "Second"), col_types = cols()) |> 
    pivot_longer(-gene, names_to = "column", values_to = "count") |>
    filter(count > 0) |>
  ggplot(aes(x = log10(count), fill = column, group = column)) +
    theme_bw() +
    geom_density(alpha = 0.3)
}

star_col_count <- function(file, smpl) {
  read_tsv(file, skip = 4, col_names = c("gene", "unstranded", "first", "second"), col_types = cols()) |> 
    filter(unstranded > 10) |> 
    summarise(unstranded = sum(unstranded), first = sum(first), second = sum(second)) |> 
    mutate(raw_sample = smpl)
}

star_strand <- function(paths, meta, suffix = ".txt", prop.limit = 0.8) {
  path <- paths$readcount
  s2n <- set_names(meta$sample, meta$raw_sample)
  m <- meta$raw_sample |> 
    map_dfr(~star_col_count(file.path(path, paste0(.x, suffix)), .x)) |> 
    mutate(sample = as.character(s2n[raw_sample])) |>
    mutate(r1 = first / unstranded, r2 = second / unstranded) |> 
    mutate(strand = if_else(r1 > prop.limit, "first", if_else(r2 > prop.limit, "second", "unstranded")))
  ms <- m |> distinct(strand)
  if (nrow(ms) == 1) {
    return(ms$strand)
  } else {
    print("Different stranding detected")
    print(m)
    stop()
  }
}

# wrapper around star_strand, returns STAR column
star_strand_column <- function(...) {
  str2col <- set_names(
    c(3, 4, 2),
    c("first","second", "unstranded")
  )
  
  strand <- star_strand(...)
  str2col[strand]
}

# Select genes with at least min.count count in at least one sample
filter_star_min_count <- function(set, min.count = 10, count.column = "count_norm"){
  set$dat <- set$dat |> 
    group_by(gene_id) |> 
    mutate(good = max(get(count.column)) >= min.count) |> 
    ungroup()
    
  set$sel <- set$dat |> 
    filter(good) |> 
    pull(gene_id) |> 
    unique()

  set  
}


filter_star_samples <- function(set, expr) {
  meta <- set$metadata |> 
    filter(eval(rlang::parse_expr(expr))) |> 
    droplevels()
  smpl_sel <- as.character(meta$sample)
  
  set$dat <- set$dat |> filter(sample %in% smpl_sel)
  set$tab <- set$tab[, smpl_sel]
  set$metadata <- meta  
  
  set
}

# returns genes with zero in at least one group
# gene_id and sum of counts per group in each group
find_zeroes <- function(set, group_var = "group") {
  set$dat |>
    left_join(set$metadata, by = c("sample", "raw_sample")) |>
    mutate(group = get(group_var)) |> 
    group_by(gene_id, group) |>
    summarise(S = sum(count)) |>
    ungroup() |>
    group_by(gene_id) |>
    mutate(m = min(S)) |>
    ungroup() |>
    filter(m == 0) |> 
    pivot_wider(id_cols = gene_id, names_from = group, values_from = S)
}


merge_star_dat <- function(dat, meta, genes, columns = c("group", "replicate")) {
  dat |> 
    left_join(meta, by = "sample") |> 
    left_join(genes, by = "gene_id") |> 
    mutate(gene_symbol = na_if(gene_symbol, "")) |> 
    select(gene_id, gene_symbol, description, sample, all_of(columns), count, count_norm, rpkm)
}


plot_qualities <- function(qc) {
  qc |> 
    ggplot(aes(x = Base, y = Mean, group = root, colour = pair)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(alpha = 0.3, size = 0.6) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = "Base", y = "Mean quality score")
}

plot_cluster_qualities <- function(qc, text.size = 10) {
  tab <- qc |> 
    pivot_wider(id_cols = c(pair, Base), names_from = raw_sample, values_from = Mean) |> 
    select(-c(pair, Base)) |> 
    as.matrix()
  
  hc <- t(tab) |> 
    dist() |> 
    hclust()
  
  dendr <- ggdendro::dendro_data(hc)
  seg <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text.size),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(size = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
  
}


plot_map_qual <- function(qc, slog, base = 20) {
  qcf <- qc |> 
    filter(Base == base) |> 
    group_by(raw_sample) |> 
    summarise(qual = mean(Mean))
  slogf <- slog |> 
    filter(key == "Uniquely mapped reads %") |> 
    mutate(mapped = str_remove(value, "%") |> as.numeric())
  qcf |> 
    left_join(slogf, by = "raw_sample") |> 
  ggplot(aes(x = mapped, y = qual)) +
    theme_bw() +
    geom_point() +
    labs(y = paste("Mean read quality at", base), x = "Uniquely mapped reads")
}


denoise_star <- function(set, id_col = "gene_id") {
  tab_nr <- noisyr::noisyr(approach.for.similarity.calculation = "counts", expression.matrix = set$tab)
  dat <- tab_nr |> 
    as_tibble(rownames = "id") |> 
    pivot_longer(-id, names_to = "sample", values_to = "count") |> 
    rename(!!id_col := id)
  list(
    dat = dat,
    tab = tab_nr,
    metadata = set$metadata
  )
}





despike <- function(set) {
  set$sel <- set$sel |> 
    str_subset("SPIKE", negate = TRUE)
  
  set
}


plot_mapped_count <- function(set) {
  sm <- colSums(set$tab)
  cnt <- tibble(
    sample = names(sm),
    value = sm,
    key = "Reads counted in genes"
  )
  counts <- set$star_log |> 
    filter(key %in% c("Number of input reads", "Uniquely mapped reads number")) |> 
    mutate(key = recode(key, "Number of input reads" = "Input reads", "Uniquely mapped reads number" = "Uniquely mapped reads")) |> 
    select(sample, value, key) |> 
    mutate(value = as.numeric(value)) |> 
    bind_rows(cnt) |> 
    mutate(value = value / 1e6) |> 
    mutate(key = fct_relevel(key, c("Input reads", "Uniquely mapped reads")))
  
  sample_input <- counts |> 
    filter(key == "Input reads") |>
    rename(input = value) |> 
    arrange(input) |> 
    mutate(sample = sample |> as_factor()) |> 
    select(-key)
  
  perc <- counts |> 
    left_join(sample_input, by = "sample") |> 
    mutate(value = 100 * value / input) |> 
    select(-input)
  
  rank_by <- function(w, sel) {
    r <- w |> 
      filter(key == sel) |> 
      mutate(rank = rank(value)) |> 
      select(sample, rank)
    w |> 
      left_join(r, by = "sample")
  }
  
  dat <- bind_rows(
      counts |>
        rank_by("Input reads") |> 
        add_column(what = "Count"),
      perc |>
        rank_by("Uniquely mapped reads") |> 
        add_column(what = "Percentage") |> 
        mutate(rank = rank + 1000)
    )
  
  labs <- dat |> 
    select(sample, rank) |> 
    distinct() |> 
    mutate(rank = as.character(rank))
  
  get_labels <- function(rnk) {
    tibble(rank = rnk) |> 
      left_join(labs, by = "rank") |> 
      pull(sample)
  }
  
  dat |> 
    ggplot(aes(x = as_factor(rank), y = value, colour = key)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    geom_segment(aes(xend = as_factor(rank), yend = 0), colour = "grey90") +
    geom_point() +
    scale_colour_manual(values = okabe_ito_palette) +
    scale_x_discrete(labels = get_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~ what, scales = "free") +
    labs(x = NULL, y = "Read count (millions)", colour = "Legend")
}


plot_chrom_proportion <- function(ids) {
  ids |>
    # filter(chr %in% CHROMOSOMES) |>
    mutate(frac = count / (length )) |> 
    # mutate(chr = factor(chr, levels = CHROMOSOMES)) |> 
  ggplot(aes(x = chr, y = frac)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_col() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Chromosome", y = "Counts per chromosome length")
}

plot_ribo_fraction <- function(ids, ribochr = "BK000964.3") {
  ids |> 
    group_by(sample) |> 
    summarise(tot_count = sum(count), ribo_count = sum(count[chr == "BK000964.3"])) |> 
    mutate(ribo_prop = ribo_count / tot_count) |> 
  ggplot(aes(y = ribo_prop, x = sample)) +
    theme_bw() +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = "Fraction of counts in ribosome") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}


select_strong_genes <- function(set, limit = 100) {
  set$dat |> 
    select(gene_id, sample, count) |> 
    left_join(set$metadata, by = "sample") |> 
    group_by(gene_id, group) |> 
    summarise(min_group_count = min(count)) |> 
    ungroup() |> 
    group_by(gene_id) |> 
    summarise(best_count = max(min_group_count)) |> 
    filter(best_count > limit) |> 
    pull(gene_id) |> 
    unique()
}


save_count_data <- function(set, file, what = "count_norm") {
  set$dat |> 
    mutate(val = get(what)) |> 
    left_join(set$genes, by = "gene_id") |> 
    pivot_wider(id_cols = c(gene_id, gene_symbol), names_from = sample, values_from = val) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_tsv(file)
}

dat2mat <- function(dat, what = "count", id_col = "gene_id") {
  dat |> 
    pivot_wider(id_cols = !!id_col, names_from = sample, values_from = !!what) |> 
    column_to_rownames(all_of(id_col)) |> 
    as.matrix()
}
