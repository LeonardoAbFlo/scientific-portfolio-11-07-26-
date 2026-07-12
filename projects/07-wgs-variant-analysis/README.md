# Whole-genome sequencing and variant analysis

Three ordered workflow collections document short-read germline calling, paired
somatic calling, and long-read bacterial variant benchmarking.

## Layout

- `workflows/short-read-germline/`: read QC and trimming, BWA-MEM alignment,
  duplicate marking, BQSR, GATK calling/recalibration, annotation, and summaries.
- `workflows/somatic-paired/`: Mutect2 calling, contamination estimation,
  FilterMutectCalls, and final quality filters.
- `workflows/long-read-benchmarking/`: depth analysis, bcftools and Clair3 calls,
  replicate filtering, annotation, gene summaries, and benchmarking.

The scripts are preserved as research workflow records and have been renamed and
reordered for readability. Several contain study-machine path defaults; set the
documented root/reference variables before use. They are not intended to be run
unchanged on private or clinical data.

