#!/usr/bin/python3

import pandas as pd

# Load annotated variants
annotated_df = pd.read_csv(
    "mtb_resistance_variants.tsv",
    sep='\t', header=None
)

# BED has columns: [chrom, start, end, ref, alt, sample_GT..., gene_chrom, gene_start, gene_end, gene_name...]
gene_col = annotated_df.columns[-1]  # Adjust if needed

# Extract gene names
annotated_df['Gene'] = annotated_df[gene_col]

# Prepare DataFrame: rows=genes, cols=barcodes, values=variant counts
samples = [col for col in annotated_df.columns if "barcode" in str(annotated_df[col][0])]
gene_list = annotated_df['Gene'].unique()

# Initialize empty DataFrame
summary_df = pd.DataFrame(0, index=gene_list, columns=samples)

# Count variants per gene and sample (barcode)
for idx, row in annotated_df.iterrows():
    gene = row['Gene']
    for sample_col in samples:
        barcode_gt = row[sample_col]
        barcode, gt = barcode_gt.split(':')
        if gt in ['1', '0/1', '1/1']:  # Variant present
            summary_df.at[gene, sample_col] += 1

summary_df.reset_index(inplace=True)
summary_df.rename(columns={'index': 'Gene'}, inplace=True)

# Save summary table
summary_df.to_csv("mtb_variant_summary.tsv", sep='\t', index=False)

print("Variant summary saved to: mtb_variant_summary.tsv")