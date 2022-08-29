wgcna_prepare <- function(set) {

  # WGCNA prefers genes in columns
  tab <- t(dat2mat(set$dat, what = "rlog"))

  # select good genes
  gsg <- WGCNA::goodSamplesGenes(tab)
  tab[gsg$goodSamples, gsg$goodGenes]
}
 
wgcna_thresholds <- function(tab, powers = c(c(1:10), seq(from = 12, to = 20, by = 2))) { 
  
  sft <- WGCNA::pickSoftThreshold(tab, powerVector = powers, verbose = 5)
  
  g1 <- sft$fitIndices |> 
    ggplot(aes(x = Power, y = -sign(slope) * SFT.R.sq, label = Power)) +
    theme_bw() +
    geom_point() +
    labs(x = "Power", y = "Scale-free topology model fit, signed R^2", title = "Scale independence")
  
  g2 <- sft$fitIndices |> 
    ggplot(aes(x = Power, y = mean.k., label = Power)) +
    theme_bw() +
    geom_point() +
    labs(x = "Power", y = "Mean connectivity", title = "Mean connectivity")
  
  list(
    sft = sft,
    plt = cowplot::plot_grid(g1, g2, nrow = 1, align = "h")
  )
}


wgcna_net <- function(tab, power, max_block_size = NULL) {
  if(is.null(max_block_size))
    max_block_size <- ncol(tab)
  
  suppressPackageStartupMessages(library(WGCNA)) # needed to mask stats::cor
  WGCNA::blockwiseModules(tab,
                  power = power,
                  maxBlockSize = max_block_size,
                  minModuleSize = 30,
                  TOMType = "signed",
                  reassignThreshold = 0,
                  mergeCutHeight = 0.25,
                  numericLabels = TRUE,
                  pamRespectsDendro = FALSE,
                  verbose = 3
  )
}



wgcna_colour_enrichment <- function(net, fterms, gene2name, fdr_limit = 0.05) {
  tb <- tibble(colour = net$colors, gene_id = names(net$colors)) |>
    filter(colour > 0)
  colours <- tb$colour |> unique() |> sort()
  all_genes <- tb$gene_id |> unique()
  ontologies <- names(fterms)
  
  map_dfr(ontologies, function(ont) {
    map_dfr(colours, function(cl) {
      cat(str_glue("Enrichment for {ont} in colour {cl}\n\n"))
      sel_genes <- tb |> 
        filter(colour == cl) |> 
        pull(gene_id)
      f <- fenr::functional_enrichment(all_genes, sel_genes, fterms[[ont]], gene2name)
      if(!is.null(f)) {
        f |> 
          filter(p_adjust < fdr_limit) |> 
          add_column(ontology = ont, colour = cl)
      } else {
        NULL
      }
    })
  })
}


wgcna_network <- function(tab, net, gene2name) {
  gene_ids <- colnames(tab)
  gene_symbols <- tibble(gene_id = gene_ids) |> 
    left_join(gene2name, by = "gene_id") |> 
    pull(gene_symbol)
  # adjacency matrix
  adj_mat <- WGCNA::adjacency(tab)
  # Useful format
  edg <- WGCNA::exportNetworkToCytoscape(adj_mat, threshold  = 0.8, nodeNames = gene_ids, altNodeNames = gene_symbols)
  clrs <- tibble(colour = net$colors, gene_id = names(net$colors)) |>
    filter(colour > 0)
  edg$edgeData |> 
    as_tibble() |> 
    left_join(clrs |> rename(from_colour = colour), by = c("fromNode" = "gene_id")) |>
    left_join(clrs |> rename(to_colour = colour), by = c("toNode" = "gene_id")) |>
    filter(from_colour == to_colour) |> 
    rename(colour = from_colour) |> 
    select(-to_colour)
}