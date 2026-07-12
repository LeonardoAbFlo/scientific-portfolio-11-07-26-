#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2, 3)) {
  stop("Usage: 11_build_sample_report.r CONTIG_REPORT OUTPUT_PREFIX [READ_QC_TSV]")
}

contigs <- read_tsv(args[[1]], show_col_types = FALSE)
output_prefix <- args[[2]]
read_qc_path <- if (length(args) == 3) args[[3]] else NA_character_

required_columns <- c(
  "classification", "ST", "K_locus", "O_locus", "completeness",
  "contamination", "molecule_type", "Best_Hit_ARO", "contig_id"
)
for (column in required_columns) {
  if (!column %in% names(contigs)) contigs[[column]] <- NA_character_
}

first_valid <- function(x) {
  values <- x[!is.na(x) & as.character(x) != ""]
  if (length(values)) values[[1]] else NA
}

sample_report <- contigs %>%
  group_by(sample_id) %>%
  summarise(
    species = first_valid(classification),
    sequence_type = first_valid(ST),
    k_locus = first_valid(K_locus),
    o_locus = first_valid(O_locus),
    completeness = suppressWarnings(as.numeric(first_valid(completeness))),
    contamination = suppressWarnings(as.numeric(first_valid(contamination))),
    total_contigs = n_distinct(contig_id),
    total_plasmids = n_distinct(contig_id[molecule_type == "plasmid"]),
    plasmids_with_amr = n_distinct(contig_id[molecule_type == "plasmid" & !is.na(Best_Hit_ARO) & Best_Hit_ARO != ""]),
    amr_genes = paste(sort(unique(na.omit(Best_Hit_ARO[Best_Hit_ARO != ""]))), collapse = ","),
    .groups = "drop"
  )

if (!is.na(read_qc_path) && file.exists(read_qc_path)) {
  sample_report <- sample_report %>%
    left_join(read_tsv(read_qc_path, show_col_types = FALSE), by = "sample_id")
}

amr_matrix <- contigs %>%
  filter(!is.na(Best_Hit_ARO), Best_Hit_ARO != "") %>%
  distinct(sample_id, Best_Hit_ARO) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = Best_Hit_ARO, values_from = present, values_fill = 0)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
write_tsv(sample_report, paste0(output_prefix, "_samples.tsv"), na = "")
write_tsv(amr_matrix, paste0(output_prefix, "_amr_matrix.tsv"), na = "")
message("Wrote sample summary and AMR matrix with prefix ", output_prefix)
