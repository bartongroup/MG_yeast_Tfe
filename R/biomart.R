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
  ), mart = mart, useCache = use_cache) |> 
  dplyr::rename(
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    gene_id = ensembl_gene_id,
    gene_symbol = external_gene_name,
    ncbi_id = entrezgene_id,
    gc_content = percentage_gene_gc_content
  ) |>
    dplyr::mutate(
      description = str_remove(description, "\\s\\[.*\\]"),
      gene_symbol = na_if(gene_symbol, "")
    ) |> 
    tibble::as_tibble() |> 
    dplyr::mutate(gene_symbol = if_else(is.na(gene_symbol), gene_id, gene_symbol))
}

biomart_fetch_proteins <- function(mart, gene_ids, use_cache = TRUE) {
  getBM(attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "ensembl_peptide_id"
  ), mart = mart, filters = "ensembl_gene_id", values = gene_ids, useCache = use_cache) |> 
    dplyr::rename(
      gene_id = ensembl_gene_id,
      gene_symbol = external_gene_name,
      protein_id = ensembl_peptide_id    
    ) |> 
    tibble::as_tibble() |> 
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
  coords <- exons |>
    group_split(ensembl_gene_id) |> 
    map_dfr(function(w) {
      IRanges::IRanges(w$exon_chrom_start, w$exon_chrom_end) |>
        IRanges::reduce() |>
        as.data.frame() |> 
        as_tibble() |> 
        mutate(ensembl_gene_id = w$ensembl_gene_id[1], chr = w$chromosome_name[1])
    }) |> 
    rename(gene_id = ensembl_gene_id, length = width)
  coords |> 
    group_by(gene_id, chr) |> 
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
  ) |> 
    dplyr::rename(gene_id = ensembl_gene_id, term_id = !!sym(id)) |> 
    dplyr::filter(term_id != "") |> 
    tibble::as_tibble()
}

# Get all GO-term descriptions

bm_fetch_go_descriptions <- function(mart, use_cache = TRUE) {
  # filtering on GO-terms does not work properly, so I have to fetch all terms
  getBM(
    attributes = c("go_id", "name_1006", "namespace_1003"),
    mart = mart, useCache = use_cache) |> 
    dplyr::rename(
      term_id = go_id,
      term_name = name_1006,
      #term_description = definition_1006,
      term_domain = namespace_1003
    ) |> 
    dplyr::filter(term_id != "") |> 
    tibble::as_tibble()
}


# Fetch GO-gene and GO descriptions from biomaRt

bm_fetch_go <- function(mart, gene_ids, slim = FALSE) {
  gene2go <- bm_fetch_go_genes(mart, gene_ids, slim)
  goterms <- bm_fetch_go_descriptions(mart)
  terms <- gene2go$term_id |> unique()
  go2gene <- map(terms, function(trm) gene2go[gene2go$term_id == trm, ]$gene_id) |> 
    set_names(terms)
  
  list(
    term2gene = go2gene,
    gene2term = gene2go,
    terms = goterms
  )
}






