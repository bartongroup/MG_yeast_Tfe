run_stringdb <- function(genes, species, score_threshold=400, version="11") {
  if(length(genes) < 3) return(NULL)
  
  if(!dir.exists("string_db")) dir.create("string_db")
  string_db <- STRINGdb$new(version="11", species=species, score_threshold=score_threshold, input_directory="string_db")
  mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  cl <- string_db$get_clusters(mapped)
  en <- map(cl, ~string_db$get_enrichment(.x))
  cl_size <- cl |> map(~length(.x)) |> unlist() 
  n_good <- sum(cl_size > 2)
  if(n_good == 0) return(NULL)
  
  cl_clip <- map(1:n_good, ~cl[[.x]])
  en <- map(cl_clip, ~string_db$get_enrichment(.x))
  list(
    string_db = string_db,
    mapped = mapped,
    clusters = cl_clip,
    enrichment = en
  )
}


plot_stringdb_clusters <- function(sdb, file, min_cluster_size = 2, size = 4000) {
  if(is.null(sdb)) return(NULL)
  
  hits <- map_dfr(sdb$clusters, function(x) {tibble(n = length(x), cl = x)}) |>
    filter(n >= min_cluster_size) |>
    pull(cl) |>
    unique()
  png(file, width = size, height = size, res = 300)
  sdb$string_db$plot_network(hits)
  dev.off()
}


proteins2geneid <- function(prots, proteins) {
  p2g <- set_names(proteins$gene_id, proteins$protein_id)
  p2g[str_remove(prots, "^\\d+\\.")] |> as.character()
}

proteins2genename <- function(prots, proteins) {
  p2g <- set_names(proteins$gene_symbol, proteins$protein_id)
  p2g[str_remove(prots, "^\\d+\\.")] |> as.character() |> na.omit()
}

stringdb_gene_clusters <- function(sdb, proteins) {
  if(is.null(sdb)) return(NULL)
  
  genes <- map_chr(sdb$clusters, function(cl) {
    proteins2genename(cl, proteins) |>  str_c(collapse = ", ")
  })
  tibble(
    cluster = seq_along(genes),
    genes = genes
  )
}
