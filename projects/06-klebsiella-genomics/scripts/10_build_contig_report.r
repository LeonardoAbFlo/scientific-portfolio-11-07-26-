#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 10_build_contig_report.r TAXONOMY_TSV PROFILE_DIR OUTPUT_TSV")
}

taxonomy_path <- args[[1]]
profile_dir <- args[[2]]
output_path <- args[[3]]

read_character_tsv <- function(path) {
  read_tsv(path, col_types = cols(.default = col_character()), show_col_types = FALSE)
}

taxonomy <- read_character_tsv(taxonomy_path)
stopifnot(all(c("sample_id", "contig_id") %in% names(taxonomy)))

rgi_paths <- list.files(
  file.path(profile_dir, "rgi"), pattern = "\\.txt$", recursive = TRUE,
  full.names = TRUE
)
rgi <- map_dfr(rgi_paths, function(path) {
  sample_id <- basename(dirname(path))
  read_character_tsv(path) %>%
    mutate(
      sample_id = sample_id,
      contig_id = str_replace(.data$Contig, "^(scaffold_[0-9]+|contig_[0-9]+).*", "\\1")
    )
})

checkm_paths <- list.files(
  file.path(profile_dir, "checkm2"), pattern = "quality_report\\.tsv$",
  recursive = TRUE, full.names = TRUE
)
checkm <- map_dfr(checkm_paths, function(path) {
  sample_id <- basename(dirname(path))
  table <- read_character_tsv(path)
  name_column <- if ("Name" %in% names(table)) table$Name else if ("Sample" %in% names(table)) table$Sample else NA_character_
  tibble(
    sample_id = sample_id,
    checkm_name = name_column,
    completeness = table$Completeness,
    contamination = table$Contamination
  )
}) %>% distinct(sample_id, .keep_all = TRUE)

kleborate_paths <- list.files(
  file.path(profile_dir, "kleborate"), pattern = "\\.tsv$", full.names = TRUE
)
kleborate <- map_dfr(kleborate_paths, function(path) {
  read_character_tsv(path) %>% mutate(sample_id = tools::file_path_sans_ext(basename(path)))
})

report <- taxonomy
if (nrow(rgi) > 0) {
  report <- rgi %>% left_join(report, by = c("sample_id", "contig_id"), relationship = "many-to-many")
}
if (nrow(checkm) > 0) report <- report %>% left_join(checkm, by = "sample_id")
if (nrow(kleborate) > 0) report <- report %>% left_join(kleborate, by = "sample_id")

preferred <- c(
  "sample_id", "contig_id", "classification", "molecule_type",
  "Best_Hit_ARO", "Drug Class", "Resistance Mechanism", "AMR Gene Family",
  "completeness", "contamination", "K_locus", "O_locus"
)
report <- report %>%
  distinct() %>%
  select(any_of(preferred), everything()) %>%
  arrange(sample_id, contig_id)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_tsv(report, output_path, na = "")
message("Wrote ", nrow(report), " rows to ", output_path)
