# Nanopore mosquito mitochondrial methylome workflow

This project profiles mitochondrial 6mA from Oxford Nanopore data. It is a
portable rewrite of the original Methyloma scripts: all paths, samples, resource
limits, models, and thresholds are now explicit inputs.

## Stages

1. Dorado basecalling with a modified-base-aware model
2. Barcode demultiplexing
3. FASTQ extraction while retaining alignment tags in read comments
4. Read mapping to sample-specific, reoriented mitochondrial references
5. Modkit pileup
6. Modkit motif discovery
7. Coverage-aware site and sample summaries

## Quick start

```bash
cp config/workflow.env.example config/workflow.env
cp config/samples.example.tsv config/samples.tsv
# Edit both files, then run:
bash run_methylome_pipeline.sh --config config/workflow.env --samples config/samples.tsv
```

The sample sheet is tab-separated:

```text
sample_id\tspecies\tbarcode\treference_fasta
F1_83\tCulex quinquefasciatus\tbarcode83\t/path/to/F1_83_reoriented.fasta
```

External software (`dorado`, `samtools`, `minimap2`, and `modkit`) must be on
`PATH`. Environment management is intentionally decoupled from the workflow.

## Interpretation

The final table distinguishes missing input, no detected motif, weak evidence,
and candidate sites meeting both coverage and modification-fraction thresholds.
Thresholds are configurable and should be justified in the associated analysis.
