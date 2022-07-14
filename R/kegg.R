# KEGG gene-pathway information
# KEGG uses NCBI identifiers internally, they need to be converted into Ensembl
# However, yeast IDs are not NCBI, just standard yeast IDs

get_kegg <- function(species, bm_genes, convert_ncbi = TRUE) {
  bm <- bm_genes %>% select(gene_id, ncbi_id) %>% drop_na() %>% distinct()
  lst <- KEGGREST::keggList("pathway", species)
  terms <- tibble(
    term_id = names(lst) %>% str_remove("path:"),
    term_name = lst
  )
  pb <- progress::progress_bar$new(total = nrow(terms))
  term2gene <- map(terms$term_id, function(path_id) {
    pw <- KEGGREST::keggGet(path_id)
    pb$tick()
    if (!is.null(pw[[1]]$GENE)) {
      # KEGG list of genes is a vector with alternate NCBI/other ID and gene description
      gns <-  pw[[1]]$GENE
      ids <- gns[seq(1, length(gns) - 1, 2)]
      if (convert_ncbi) {
        ids <- bm %>%
          filter(ncbi_id %in% ids) %>%
          pull(gene_id) %>%
          unique()  # convert NCBI to Ensembl, warning: not one-to-one!
      }
      ids
    }
  }) %>% 
    set_names(terms$term_id)
  
  gene2term <- map_dfr(terms$term_id, function(tid) {
    tibble(
      gene_id = term2gene[[tid]],
      term_id = tid
    )
  })
  
  list(
    terms = terms,
    term2gene = term2gene,
    gene2term = gene2term
  )
}