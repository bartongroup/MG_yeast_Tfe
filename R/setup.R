EXPERIMENTS <- tribble(
  ~NAME, ~PATH, ~SAMPLE_FILE, ~ENSEMBL_DATASET, ~ENSEMBL_VERSION, ~EXAMPLE, ~STAR_COLUMN,
  "tfe", "rna_seq", "rna_seq/config/samples.txt", "scerevisiae_gene_ensembl", "106", "165-1", 2
)

PATH <-  "rna_seq"
SAMPLE_FILE <- "rna_seq/config/samples.txt"
ENSEMBL_DATASET <- "scerevisiae_gene_ensembl"
ENSEMBL_VERSION <- "106"
KEGG_SPECIES <- "sce"
GO_SPECIES <- "sgd"

EXAMPLE <- "WT_0_1"

LOGFC_LIMIT <- 1
FDR_LIMIT <- 0.01

CONTRASTS_0 <- c("Tfe2_0-WT_0", "Tfe1_0-WT_0", "Tfe1_0-Tfe2_0")
CONTRAST_SELECTION <- c("Tfe1_0-WT_0", "Tfe1_30-WT_30", "Tfe1_60-WT_60", "Tfe2_0-WT_0", "Tfe2_30-WT_30", "Tfe2_60-WT_60")


SAMPLE_RENAME <- tribble(
  ~raw_sample, ~strain, ~time, ~replicate,
  "WT-60-1", "WT", "60", "1",
  "WT-60-2", "Tfe2", "60", "2",   # 
  "WT-60-3", "Tfe2", "60", "3",   #
  "167-60-1", "Tfe2", "60", "1",     #
  "167-60-2", "Tfe1", "60", "2",
  "167-60-3", "Tfe1", "60", "3",
  "165-60-1", "Tfe1", "60", "1",     #
  "165-60-2", "WT", "60", "2",     #
  "165-60-3", "WT", "60", "3"   #
)


make_metadata <- function(name, sample_file, sample_rename) {
  if (name == "tfe") {
    meta <- read_tsv(sample_file, col_names = "raw_sample", show_col_types = FALSE) %>%
      separate(raw_sample, c("strain_id", "time", "replicate"), sep = "-", fill = "right", remove = FALSE) %>% 
      mutate(
        mis = is.na(replicate),
        replicate = if_else(mis, time, replicate),
        time = if_else(mis, "0", time),
        strain = recode(strain_id, "165" = "Tfe2", "167" = "Tfe1")
      ) %>% 
      left_join(sample_rename, by = "raw_sample", suffix = c("", "_replace")) %>% 
      mutate(scramble = !is.na(strain_replace) & (strain != strain_replace | replicate != replicate_replace)) %>% 
      unite(raw_group, c(strain_id, time), sep = "_", remove = FALSE) %>% 
      mutate(
        strain = if_else(is.na(strain_replace), strain, strain_replace),
        time = if_else(is.na(time_replace), time, time_replace),
        replicate = if_else(is.na(replicate_replace), replicate, replicate_replace)
      ) %>% 
      select(-ends_with("replace")) %>% 
      unite(sample, c(strain, time, replicate), sep = "_", remove = FALSE) %>% 
      unite(group, c(strain, time), sep = "_", remove = FALSE) %>% 
      mutate(
        strain_id = as_factor(strain_id) %>% fct_relevel(c("WT", "167")),
        strain = as_factor(strain) %>% fct_relevel(c("WT", "Tfe1"))
      ) %>% 
      arrange(strain, time, replicate) %>% 
      mutate(
        across(c(sample, time, replicate, group), as_factor)
      ) %>% 
      select(raw_sample, sample, strain_id, strain, raw_group, group, time, replicate, scramble)
  }
  meta
}

make_dirs <- function(top_dir) {
  list(
    starmap = file.path(top_dir, "starmap"),
    fscreen = file.path(top_dir, "fscreen"),
    readcount = file.path(top_dir, "readcount"),
    chrcount = file.path(top_dir, "chrcount"),
    bedgraph = file.path(top_dir, "bedgraph"),
    qc = file.path(top_dir, "qc")
  )
}

# Gene names in the W303 genome have a suffix "_W303" or "mRNA". Need to remove
# these before processing.
fix_gene_names <- function(v, suffixes = c("_W303", "_mRNA")) {
  expr <- str_c(suffixes, collapse = "|")
  v %>% 
    str_remove(expr)
}
