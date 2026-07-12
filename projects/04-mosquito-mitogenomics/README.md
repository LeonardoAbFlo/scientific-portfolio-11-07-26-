# Mosquito mitogenomics workflow

A modular long-read workflow for mosquito mitochondrial genome reconstruction and
comparative analysis.

## Modules

- `workflow/01-read-processing-and-assembly/`: basecalling, demultiplexing,
  FASTQ/QC generation, Flye assembly selection, mapping, and polishing.
- `workflow/02-annotation-and-comparative-genomics/`: MitoFinder/MITOS2
  annotation, manual-correction handoff, genome statistics, codon usage,
  amino-acid composition, Circos, and synteny tables.
- `workflow/03-phylogenetics/`: reference annotation, dataset preparation,
  concatenation of 13 mitochondrial protein-coding genes, tree inference, and
  rooting.

The modified-base module was separated and fully rewritten as the neighboring
`05-mosquito-methylome` project.

These scripts retain their stage numbers and environment hooks because they
encode an end-to-end research workflow. Configure the workspace and reference
paths before execution; large databases and sequence data are not included.

