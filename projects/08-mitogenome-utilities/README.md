# Mitogenome utilities

Supporting scripts developed around mosquito mitochondrial assembly and
annotation:

- `trna-annotation/`: ARWEN execution and tRNA extraction
- `mitofinder/`: final trnI- and rrnS-based MitoFinder/Medaka workflows
- `reorientation/`: final control-region and trnI sequence-reorientation methods
- `dnaapler/`: canonical DNAapler-based circular-sequence reorientation
- `assembly-graph/`: Bandage graph extraction

Canonical scripts selected from historical version families are explicitly
labeled:

- `dnaapler/reorient_mitogenomes_with_dnaapler_final_version.sh`
- `mitofinder/reorient_by_trni_and_reannotate_final_version.sh`
- `mitofinder/reorient_by_rrns_and_reannotate_final_version.sh`
- `reorientation/reorient_by_trni_final_version.sh`

Only the final implementation of each method is retained. Reference databases
and sequence files are intentionally absent.
