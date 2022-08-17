SPECIES <- list(
  human = list(
    go = "goa_human",
    re = "Homo sapiens",
    kg = "hsa"
  ),
  yeast = list(
    go = "sgd",
    re =  "Saccharomyces cerevisiae",
    kg = "sce"
  )
)


get_functional_terms <- function(geneids, species, all_features) {
  sp <- SPECIES[[species]]
  cat("Loading GO term data\n")
  go <- fenr::fetch_go_from_go(sp$go)
  cat("Loading Reactome data\n")
  re <- fenr::fetch_reactome(sp$re)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(sp$kg)
  cat("Loading BioPlanet data\n")
  bp <- fenr::fetch_bp()
  
  go$mapping <- go$mapping |> 
    rename(gene_id = gene_synonym)
  bp$mapping <- bp$mapping |>
    inner_join(select(geneids, gene_id, gene_symbol), by = "gene_symbol")

  terms <- list(
    go = go,
    re = re,
    kg = kg,
    bp = bp
  )
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features)
  }) |> 
    set_names(ontologies)
}