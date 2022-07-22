sh_plot_volcano <- function(d, alpha = 0.05, title = NULL) {
  sres <- d %>% filter(FDR <= alpha)
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
  brkn <- c(-brkn, brkn[2:length(brks)]) %>% sort()
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


sh_plot_genes <- function(dat, meta, genes, scale, max_points = 100) {
  dat <- dat %>% 
    filter(gene_id %in% genes) 
  n <- nrow(dat)
  if (n == 0) return(NULL)
  
  dat <- dat %>% 
    left_join(meta, by = c("sample", "group", "replicate"))
  
  if (scale == "log") {
    dat$rpkm <- log10(dat$rpkm)
  }
  d <- dat %>%
    group_by(sample) %>% 
    summarise(
      m = mean(rpkm, na.rm = TRUE),
      se = sd(rpkm, na.rm = TRUE) / sqrt(sum(!is.na(rpkm))),
      group = first(group),
      replicate = first(replicate) %>% as_factor(),
      time = first(time),
      strain = first(strain)
    ) %>% 
    replace_na(list(se = 0)) %>% 
    mutate(lo = m - se, up = m + se) %>% 
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
  if (scale == "log") {
    g <- g + labs(y = "Log RPKM")
  } else {
    g <- g +
      labs(y = "RPKM") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06)), limits = c(0, NA))
  }
  if (n > 1) {
    g <- g + geom_errorbar(data = d, aes(x = time, ymin = lo, ymax = up, group = replicate), colour = "grey70", alpha = 0.5, position = pd, width = 0.1)
  }
  g
}




# create a few structure for fast selection in enrichment
sh_prepare_for_enrichment <- function(data, universes = c("go", "re", "kg")) {
  genes_all <- data$genes$gene_id %>% unique()
  
  map(universes, function(u) {
    term_data <- data$terms[[u]]
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
        filter(gene_id %in% genes_all)
      
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



# New version, needs preparation (see above)
sh_functional_enrichment <- function(genes_all, genes_sel, term_data, gene2name = NULL,
                                     min_count = 2, sig_limit = 0.05) {
  
  # all terms present in the selection
  terms <- term_data$gene_term %>% 
    filter(gene_id %in% genes_sel) %>% 
    pull(term_id) %>% 
    unique()
  
  # number of selected genes
  Nsel <- length(genes_sel)
  # size of the universe
  Nuni <- length(genes_all)
  
  # empty line for missing terms
  na_term <- term_data$term_info[[1]] %>% mutate_all(~NA_character_)
  
  res <- map_dfr(terms, function(term) {
    info <- term_data$term_info[[term]]
    # returns NAs if no term found
    if (nrow(info) == 0) info <- na_term %>% mutate(term_id = term)
    
    # all genes with the term
    tgenes <- term_data$term2gene[[term]]
    
    # genes from selection with the term
    tgenes_sel <- intersect(tgenes, genes_sel)
    
    nuni <- length(tgenes)
    nsel <- length(tgenes_sel)
    
    expected <- nuni * Nsel / Nuni
    # fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow = 2)
    # ft <- fisher.test(fish, alternative = "greater")
    # p <- as.numeric(ft$p.value)
    
    # Hypergeometric function much fater than fisher.test
    p <- 1 - phyper(nsel - 1, nuni, Nuni - nuni, Nsel)
    
    if (!is.null(gene2name)) tgenes_sel <- gene2name[tgenes_sel] %>% unname()
    
    bind_cols(
      info,
      tibble(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        ids = paste(tgenes_sel, collapse = ", "),
        P = p
      )
    )
  }) %>% 
    mutate(P = p.adjust(P, method = "BH")) %>% 
    filter(sel >= min_count & P <= sig_limit) %>% 
    arrange(desc(enrich)) %>% 
    mutate(enrich = round(enrich, 1), expect = round(expect, 2))
  
  res
}


