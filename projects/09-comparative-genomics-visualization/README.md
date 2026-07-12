# Comparative-genomics visualization

Python analyses for mitochondrial synteny, Circos-style genome maps, relative
synonymous codon usage, and amino-acid composition.

Recommended entry points:

- `synteny/mitogenome_synteny_by_genus_final_version.py`
- `circos/comparative_circular_genome_plot.py`
- `circos/mitogenome_circos_by_genus_final_version.py`
- `codon-usage/plot_codon_usage_by_genus_final_version.py`
- `amino-acid-composition/plot_amino_acid_composition_by_genus_final_version.py`

Superseded iterations and one-sample prototypes were removed. Each analytical
task now has one canonical final implementation. Scripts selected from numbered
historical revision families carry the literal `_final_version` label.

Most recent scripts expose command-line paths; where historical defaults remain,
pass project and reference directories explicitly. Generated figures and input
sequence data are not included.
