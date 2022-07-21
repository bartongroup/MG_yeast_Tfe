
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




read_reactome_file <- function(file_name, colnames,  reactome_url = "https://reactome.org/download/current") {
  url <- str_glue("{reactome_url}/{file_name}")
  read_tsv(url, col_names = colnames, show_col_types = FALSE)
}

reactome_fetch_pathways <- function() {
  cache_file <- "cache/reactome_pathways.tsv"
  if (!file.exists(cache_file)) {
    read_reactome_file("ReactomePathways.txt", c("reactome_id", "name", "species")) %>% 
      write_tsv(cache_file)
  }
  read_tsv(cache_file, show_col_types = FALSE)
}

reactome_fetch_genes <- function(species) {
  sp <- str_replace(species, "\\s", "_") %>% tolower()
  cache_file <- str_glue("cache/ensembl2reactome_{sp}.tsv")
  if (!file.exists(cache_file)) {
    read_reactome_file("Ensembl2Reactome.txt", c("gene_id", "term_id", "url", "event", "evidence", "species_")) %>% 
      filter(species_ == species) %>% 
      select(gene_id, term_id) %>% 
      distinct() %>% 
      write_tsv(cache_file)
  }
  read_tsv(cache_file, show_col_types = FALSE)
}

fetch_reactome <- function(species, bm_genes) {
  r <- reactome_fetch_pathways()
  g2r <- reactome_fetch_genes(species) %>% 
    left_join(select(bm_genes, gene_id, gene_name), by = "gene_id") %>% 
    drop_na() %>% 
    select(gene_id, term_id)
  terms <- g2r$term_id %>%
    unique()
  
  reactometerms <- select(r, reactome_id, name) %>%
    rename(term_id = reactome_id, term_name = name) %>%
    filter(term_id %in% terms) %>%
    distinct()
  r2g <- map(terms, function(trm) g2r[g2r$term_id == trm, ]$gene_id) %>%
    set_names(terms)
  list(
    gene2term = g2r,
    term2gene = r2g,
    terms = reactometerms
  )
}



# Get ontology directly from geneontology.org

go_fetch_go_descriptions <- function(obo_file = "cache/go.obo") {
  if (!file.exists(obo_file)) {
    download.file("http://purl.obolibrary.org/obo/go.obo", obo_file)
  }
  go <- ontologyIndex::get_ontology(obo_file, extract_tags = c("everything"))
  tibble(
    term_id = go$id,
    term_name = go$name,
    term_namespace = unlist(go$namespace)
  )
}

# GO annotations from geneontology.org

GAF_COLUMNS <- c(
  "db", "db_id", "symbol", "qualifier", "go_term", "db_ref", "evidence", "from", "aspect", "db_object_name", "db_object_synonym", "db_object_type", "taxon", "date", "assigned_by", "annotation_extension", "form_id"
)

go_fetch_go_genes <- function(species, gaf_file = "cache/go_gene.gaf.gz") {
  if (!file.exists(gaf_file)) {
    download.file(str_glue("http://current.geneontology.org/annotations/{species}.gaf.gz"), gaf_file)
  }
  read_tsv(gaf_file, comment = "!", quote = "", col_names = GAF_COLUMNS, show_col_types = FALSE) %>% 
    mutate(gene_id = str_remove(db_object_synonym, "\\|.*$")) %>% 
    select(gene_id, term_id = go_term) %>% 
    distinct()
}


# Fetch GO-gene and GO descriptions from geneontolgy.org

go_fetch_go <- function(species) {
  gene2go <- go_fetch_go_genes(species)
  goterms <- go_fetch_go_descriptions()
  terms <- gene2go$term_id %>% unique()
  go2gene <- map(terms, function(trm) gene2go[gene2go$term_id == trm, ]$gene_id) %>% 
    set_names(terms)
  
  list(
    term2gene = go2gene,
    gene2term = gene2go,
    terms = goterms
  )
  
}


