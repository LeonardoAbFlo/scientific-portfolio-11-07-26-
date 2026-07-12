# Klebsiella long-read genomics workflow

This project is a reproducible rewrite of a research workflow for Oxford Nanopore
assemblies from *Klebsiella* and related Enterobacterales. It assembles reads,
profiles plasmids and antimicrobial-resistance genes, evaluates genome quality and
taxonomy, filters mixed assemblies with an auditable decision table, and builds
contig- and sample-level reports.

## Workflow

1. `scripts/01_assemble_reads.sh` runs Flye for each FASTQ file.
2. `scripts/02_split_contigs.py` writes one FASTA per contig and a manifest.
3. `scripts/03_profile_assemblies.sh` runs MOB-suite and RGI on the unfiltered assemblies.
4. `scripts/04_classify_contigs.sh` runs GTDB-Tk and CheckM2 on split contigs.
5. `scripts/05_build_taxonomy_table.py` combines contig, plasmid, and taxonomy results.
6. `scripts/06_filter_contigs.py` compares observed taxonomy with expected metadata and
   copies retained contigs to a new directory. It never deletes source data.
7. `scripts/07_merge_contigs.py` rebuilds one filtered assembly per sample.
8. `scripts/08_profile_filtered_assemblies.sh` runs MOB-suite, RGI, CheckM2, and Kleborate.
9. `scripts/09_assess_read_quality.sh` extracts NanoStat metrics from the original reads.
10. `scripts/10_build_contig_report.r` joins tool outputs into a contig-level report.
11. `scripts/11_build_sample_report.r` creates a sample summary and AMR presence matrix.

`run_klebsiella_pipeline.sh` coordinates the computational stages. Reporting is kept
separate because published studies often require manual review of tool outputs
and sample metadata before tables are frozen.

## Quick start

```bash
cp config/workflow.env.example config/workflow.env
# Edit paths and resource settings in config/workflow.env.
bash run_klebsiella_pipeline.sh --config config/workflow.env
```

Each external tool should be available on `PATH`. Conda, Apptainer, or module
activation is intentionally left to the execution environment instead of being
hard-coded in the scripts.

## Metadata contract

The filtering step accepts a tab-separated metadata file with:

```text
sample_id\texpected_species
GRC001\tKlebsiella pneumoniae
```

The taxonomy table must contain sample, contig, molecule-type, and classification
columns. Common historical names (`barcode`, `UCL_ID`, `contig_id`, `circ`) are
recognized automatically. All exclusions are written to `filter_decisions.tsv`.

## Scientific safeguards

- Source assemblies and split contigs are read-only.
- Unknown classifications are retained for review.
- Chromosome exclusions require a confident genus mismatch.
- Cross-genus plasmid filtering is opt-in because plasmid host assignment is
  uncertain and horizontal transfer is biologically plausible.
- Software versions and logs are written under the results directory.
