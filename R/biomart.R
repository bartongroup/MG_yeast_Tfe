biomart_fetch_genes <- function(mart, use_cache = TRUE) {
  getBM(attributes = c(
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "percentage_gene_gc_content",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id",
    "description"
  ), mart = mart, useCache = use_cache) %>% 
  dplyr::rename(
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    gene_id = ensembl_gene_id,
    gene_name = external_gene_name,
    ncbi_id = entrezgene_id,
    gc_content = percentage_gene_gc_content
  ) %>%
    dplyr::mutate(
      description = str_remove(description, "\\s\\[.*\\]"),
      gene_name = na_if(gene_name, "")
    ) %>% 
  tibble::as_tibble()
}

biomart_fetch_proteins <- function(mart, gene_ids, use_cache = TRUE) {
  getBM(attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_peptide_id"
  ), mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) %>% 
    dplyr::rename(
      gene_id = ensembl_gene_id,
      gene_name = external_gene_name,
      protein_id = ensembl_peptide_id    
    ) %>% 
    tibble::as_tibble() %>% 
    filter(protein_id != "")
}

# this is inspired by getGeneLengthAndGCContent function in EDASeq
biomart_fetch_gene_lengths <- function(mart, gene_ids, use_cache = TRUE) {
  exons <- getBM(attributes = c(
    "ensembl_gene_id",
    "ensembl_exon_id",
    "chromosome_name",
    "exon_chrom_start",
    "exon_chrom_end"),
    mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache
  )
  # find all exons, reduce overlap
  coords <- exons %>%
    group_split(ensembl_gene_id) %>% 
    map_dfr(function(w) {
      IRanges::IRanges(w$exon_chrom_start, w$exon_chrom_end) %>%
        IRanges::reduce() %>%
        as.data.frame %>% 
        as_tibble %>% 
        mutate(ensembl_gene_id = w$ensembl_gene_id[1], chr = w$chromosome_name[1])
    }) %>% 
    rename(gene_id = ensembl_gene_id, length = width)
  coords %>% 
    group_by(gene_id, chr) %>% 
    summarise(length = sum(length))
}


# For a given list of gene names, fetch all GO-terms

bm_fetch_go_genes <- function(mart, gene_ids, slim = FALSE, use_cache = TRUE) {
  id <- ifelse(slim, "goslim_goa_accession", "go_id")
  gene2go <- getBM(
    attributes = c("ensembl_gene_id", id),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart,
    useCache = use_cache
  ) %>% 
    dplyr::rename(gene_id = ensembl_gene_id, term_id = !!sym(id)) %>% 
    dplyr::filter(term_id != "") %>% 
    tibble::as_tibble()
}

# Get all GO-term descriptions

bm_fetch_go_descriptions <- function(mart, use_cache = TRUE) {
  # filtering on GO-terms does not work properly, so I have to fetch all terms
  getBM(
    attributes = c("go_id", "name_1006", "namespace_1003"),
    mart = mart, useCache = use_cache) %>% 
    dplyr::rename(
      term_id = go_id,
      term_name = name_1006,
      #term_description = definition_1006,
      term_domain = namespace_1003
    ) %>% 
    dplyr::filter(term_id != "") %>% 
    tibble::as_tibble()
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


# Fetch GO-gene and GO descriptions from biomaRt

bm_fetch_go <- function(mart, gene_ids, slim = FALSE) {
  gene2go <- bm_fetch_go_genes(mart, gene_ids, slim)
  goterms <- bm_fetch_go_descriptions(mart)
  terms <- gene2go$term_id %>% unique()
  go2gene <- map(terms, function(trm) gene2go[gene2go$term_id == trm, ]$gene_id) %>% 
    set_names(terms)
  
  list(
    term2gene = go2gene,
    gene2term = gene2go,
    terms = goterms
  )
}




# Reactome

reactome_fetch_genes <- function(mart, gene_ids, use_cache = TRUE) {
  gene2re <- getBM(
    attributes = c("ensembl_gene_id", "reactome"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart,
    useCache = use_cache
  ) %>% 
    dplyr::rename(gene_id = ensembl_gene_id, term_id = "reactome") %>% 
    dplyr::filter(term_id != "") %>% 
    tibble::as_tibble()
}


reactome_fetch_pathways <- function() {
  url <- "https://reactome.org/download/current/ReactomePathways.txt"
  colms <- c("reactome_id", "name", "species")
  read_tsv(url, col_names = colms, col_types = cols())
}



# Get Reactome data in the same format as GO-data

fetch_reactome <- function(mart, gene_ids) {
  r <- reactome_fetch_pathways()
  g2r <- reactome_fetch_genes(mart, gene_ids)
  terms <- g2r$term_id %>% unique()
  
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



