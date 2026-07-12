#!/usr/bin/env Rscript

# load libraries
suppressWarnings(suppressMessages(library(data.table)))

# 1) File paths
bed_file     <- "/path/to/Desktop/benchmarking/references/tbdb.bed"   # your BED file
depth_folder <- "/path/to/Desktop/benchmarking/depth_results/processed_data/depth/stop_receiving" # folder containing depth_barcode01.txt, etc.
output_file  <- "/path/to/Desktop/benchmarking/tables/pass_depth_1.csv"

# 2) List depth files, selecting only the first 8 barcodes
depth_files <- list.files(depth_folder, 
                          pattern = "^depth_barcode.*\\.txt$", 
                          full.names = TRUE)[8]  # Select only the first 8 files

# 3) Read BED file and add gene labels (assuming the 4th column is gene ID/name)
dt_bed <- fread(bed_file,
                col.names = c("chromosome","start","end","gene"),
                select = 1:4)
setkey(dt_bed, chromosome, start, end)

# 4) Read each depth file, tag with 'barcode' name, combine
all_depths <- list()
for (f in depth_files) {
  tmp <- fread(f, col.names = c("chromosome","position","depth"))
  tmp[, barcode := gsub("\\.txt$", "", basename(f))]
  all_depths[[f]] <- tmp
}
dt_depth_all <- rbindlist(all_depths, use.names=TRUE, fill=TRUE)

# 5) Prepare for interval overlap
dt_depth_all[, start := position]
dt_depth_all[, end   := position]
setkey(dt_depth_all, chromosome, start, end)

# 6) Interval join using data.table::foverlaps
merged <- foverlaps(dt_depth_all, dt_bed, 
                    by.x=c("chromosome","start","end"),
                    type="any", nomatch=NA)

# 7) Label clearly on/off-target and gene name
merged[, target := ifelse(is.na(gene), "off_target", "on_target")]

# 8) Clean up columns clearly
final_table <- merged[, .(
  chromosome = chromosome,
  position   = position,
  depth      = depth,
  barcode    = barcode,
  target     = target,
  gene       = gene
)]

# 9) Write out
fwrite(final_table, output_file, sep="\t")

# Check output
print(head(final_table))

# ---------------------------------------------------------------------
# At this point, final_table has all barcodes combined and labeled.
# You can now plot in R (ggplot2 or otherwise).
# ---------------------------------------------------------------------

# Example plotting (uncomment if you run interactively):
# library(ggplot2)
# ggplot(final_table, aes(x=position, y=depth, color=on_target)) +
#   geom_point() +
#   facet_wrap(~ barcode) +
#   theme_minimal() +
#   labs(x="Genome Position", 
#        y="Depth", 
#        title="Depth by Barcode, On vs Off-Target")