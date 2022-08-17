parse_one_fscreen <- function(file, raw_sample) {
  read_tsv(file, skip = 1, comment = "%Hit", col_types = cols()) |> 
    pivot_longer(-Genome) |> 
    mutate(raw_sample = raw_sample)
}


parse_fscreens <- function(paths, meta, suffix = "_R1_screen.txt") {
  path <- paths$fscreen
  meta$raw_sample |> 
    map_dfr(~parse_one_fscreen(file.path(path, paste0(.x, suffix)), .x)) |> 
    left_join(select(meta, raw_sample, sample), by = "raw_sample")
}

plot_fscreen_map <- function(fscreen) {
  fscreen |> 
    filter(name == "%One_hit_one_genome") |> 
    mutate(value = na_if(value, 0)) |> 
    mutate_at(c("Genome", "sample"), as_factor) |> 
    ggplot(aes(x = sample, y = Genome, fill = value)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    geom_tile() +
    scale_fill_viridis(option = "cividis", trans = "log10", limits = c(0.01, 100), labels = c(0.01,0.1,1,10,100)) +
    labs(x = NULL, y = "Genome", fill = "% Unique hits")
}


plot_fscreen_sample <- function(fscreen, smpl) {
  fscreen |>
    filter(sample == smpl & str_detect(name, "%") & !str_detect(name, "Unma")) |>
    mutate(Genome = as_factor(Genome) |> fct_reorder(value, sum)) |> 
    mutate(name = str_replace_all(name, "_", " ")) |> 
    ggplot(aes(x = Genome, y = value, fill = name)) +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette[c(1,2,3,5)]) +
    coord_flip() +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Genome", y = "Percentage", title = smpl, fill = "Quantity")
  
}

# for paired
make_roots_paired <- function(meta, suffix = "_R") {
  roots_1 <- paste0(meta$raw_sample, suffix, "1")
  roots_2 <- paste0(meta$raw_sample, suffix, "2")
  tibble(
    root = c(roots_1, roots_2),
    raw_sample = rep(meta$raw_sample, 2),
    pair = c(rep(1, nrow(meta)), rep(2, nrow(meta))) |> as_factor()
  )
}

# for single
make_roots_single <- function(meta) {
  roots <- meta$raw_sample
  tibble(
    root = roots,
    raw_sample = meta$raw_sample,
    pair = rep(1, nrow(meta)) |> as_factor()
  )
}

parse_one_qc <- function(zip_file, root) {
  if (!file.exists(zip_file)) {
    warning(paste("File", zip_file, "not found."))
    return(NULL)
  }
  dirname <- str_remove(basename(zip_file), ".zip")
  con <- unz(zip_file, file.path(dirname, "fastqc_data.txt"), open = "rb")
  s <- read_lines(con)
  close(con)
  p1 <- str_which(s, "Per base sequence quality")[1] + 1
  p2 <- str_which(s, "END_MODULE")
  p2 <- min(p2[p2 > p1]) - 1
  st <- s[p1:p2] |> str_remove("#")
  read_tsv(I(st), show_col_types = FALSE, progress = FALSE) |>
    mutate(Base = str_remove(Base, "\\-.*$") |> as.integer()) |> 
    mutate(root = root)
}

parse_qcs <- function(paths, meta, paired, suffix = "_fastqc.zip") {
  if (paired) {
    roots <-  make_roots_paired(meta)
  } else {
    roots <-  make_roots_single(meta)
  }
  path <- paths$qc
  roots$root |> 
    map_dfr(~parse_one_qc(file.path(path, paste0(.x, suffix)), .x)) |> 
    left_join(roots, by = "root")
}


plot_fscreen_slog <- function(fs, slog, genome, fs_var = "%One_hit_one_genome", slog_var = "Uniquely mapped reads %") {
  fsf <- fs |> 
    filter(name == fs_var & Genome == genome) |> 
    group_by(raw_sample) |> 
    summarise(screen = mean(value))
  slogf <- slog |> 
    filter(key == slog_var) |> 
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    rename(star = value)
  fsf |> 
    left_join(slogf, by = "raw_sample") |> 
  ggplot(aes(x = star, y = screen)) +
    theme_bw() +
    geom_point() +
    labs(x = slog_var, y = fs_var)
}




plot_flow_qc <- function(slog, qc, meta, base = 20) {
  qcf <- qc |> 
    filter(Base == base) |> 
    group_by(raw_sample) |> 
    summarise(quality = mean(Mean)) 
  slogf <- slog |> 
    filter(key == "Uniquely mapped reads %") |> 
    mutate(mapped_prop = str_remove(value, "%") |> as.numeric()) |> 
    select(raw_sample, mapped_prop)
  meta |> 
    left_join(qcf, by = "raw_sample") |> 
    left_join(slogf, by = "raw_sample") |> 
    pivot_longer(cols = c(quality, mapped_prop)) |> 
  ggplot(aes(x = value, y = flow_cell)) +
    theme_bw() +
    geom_beeswarm(groupOnX = FALSE, size = 0.5, cex = 0.5) +
    facet_wrap(~name, scales = "free_x") +
    labs(x = NULL, y = "Flow cell")
}

# read file consisting of first fastq lines, created with
# for file in *gz; do printf "$file "; gunzip -c $file | head -n 1; done > fastq_headers.txt
read_fastq_headers <- function(file) {
  read_delim(file, " ", col_names = c("file_name", "id", "attr"), col_types = "ccc") |> 
    mutate(raw_sample = str_remove(file_name, "_[12].fq.gz")) |> 
    mutate(accession = str_remove(id, "@") |> str_remove("\\.1")) |> 
    separate(attr, c("flow_cell", "lane", "tile", "x", "yr"), sep = ":") |> 
    separate(yr, c("y", "pair"), sep = "/") |> 
    group_by(raw_sample) |> 
    summarise(flow_cell = first(flow_cell) |> as_factor(), accession = first(accession) |> as.factor())
}
