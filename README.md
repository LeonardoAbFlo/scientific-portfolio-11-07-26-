# Computational biology and biomedical data science portfolio

This repository presents research software spanning genomics, molecular
simulation, machine learning, and neuroinformatics. The projects were reorganized
from working research directories into topic-focused, reviewable workflows for an
academic application. Raw data, generated images, logs, binaries, and transient
analysis files are intentionally excluded.

## Projects

| Project | Research focus | Main technologies |
|---|---|---|
| [RNA hairpin molecular dynamics](projects/01-rna-hairpin-molecular-dynamics/) | GROMACS setup, simulation, trajectory analysis, and structural summaries | Bash, Python, GROMACS |
| [GBA virtual screening](projects/02-gba-virtual-screening/) | ChEMBL curation, structure preparation, GNINA docking, and affinity modeling | Python, RDKit, scikit-learn, GNINA |
| [Stroke rehabilitation EEG](projects/03-stroke-rehabilitation-eeg/) | Leakage-safe motor-imagery decoding and electrode contribution analysis | Python, R, MNE-style signal processing, scikit-learn |
| [Mosquito mitogenomics](projects/04-mosquito-mitogenomics/) | Nanopore assembly, mitochondrial annotation, comparative genomics, and phylogenetics | Bash, Python, Flye, MitoFinder, IQ-TREE |
| [Mosquito methylome](projects/05-mosquito-methylome/) | Modified-base-aware basecalling, mitochondrial 6mA profiling, and motif discovery | Bash, Python, Dorado, minimap2, Modkit |
| [Klebsiella genomics](projects/06-klebsiella-genomics/) | Long-read assembly, plasmids, AMR, taxonomy, contamination filtering, and reporting | Bash, Python, R, GTDB-Tk, MOB-suite, RGI, Kleborate |
| [WGS and variant analysis](projects/07-wgs-variant-analysis/) | Germline, paired somatic, and long-read variant workflows | Bash, Python, R, GATK, bcftools, Clair3 |
| [Mitogenome utilities](projects/08-mitogenome-utilities/) | tRNA annotation, reorientation, assembly graphs, and MitoFinder helpers | Python, Bash, ARWEN, DNAapler |
| [Comparative-genomics visualization](projects/09-comparative-genomics-visualization/) | Synteny, Circos, codon usage, and amino-acid composition | Python, Biopython, pandas, matplotlib |

## Reproducibility notes

Many workflows depend on scientific command-line software and reference
databases that cannot reasonably be committed to Git. Each project documents its
expected tools and inputs. Results are not included as evidence unless their
underlying data can be shared ethically and legally.
