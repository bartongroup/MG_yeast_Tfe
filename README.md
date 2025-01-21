# Tfe1 and Tfe2 toxin expression in S. cerevisiae RNA-seq

Collaborators: Katharina Trunk

## Usage

### On a Linux cluster

Create and activate a `conda` environment

```
cd rna_seq
conda create --name yeast_tfe --file env.txt
conda activate yeast_tfe
```

Make sure FASTQ files are in the ./fastq subdirectory.

Run snakemake

```
./run_snake.sh
```

This will trim adapters, perform quality control, download genome files, map reads to the reference and count reads per gene.

### In RStudio

Once snakemake is finished, we suggest using RStudio. If this is done on a different machinge (I run RStudio on a laptop), some data need to be copied over (see `./get_data.sh` and `./rsync_include.txt`). Once in RStudio, start in the top project directory. The first step is to create environment using `renv`:

```
install.packages("renv")
renv::restore()
```

This will install all necessary packages. Run the `targets` pipeline.

```
targets::tar_make()
```

This will carry out all the calculations, create figures (some as [targets](https://books.ropensci.org/targets/), some in `./fig` directory) and output TSV files in directory `./tab`.


