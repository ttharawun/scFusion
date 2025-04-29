# Tint Created for Mapping Barcodes with Reads

import pysam
import sys

# Take command line arguments
bam_file = sys.argv[1]
chimeric_sam = sys.argv[2]
output_sam = sys.argv[3]

# Step 1: Build dictionary of {read_name: CB_tag}
print(f"Indexing BAM file {bam_file} for read → CB mapping...")
bam = pysam.AlignmentFile(bam_file, "rb")
cb_dict = {}
for read in bam:
    cb = read.get_tag("CB") if read.has_tag("CB") else None
    if cb:
        cb_dict[read.query_name] = cb
bam.close()

# Step 2: Annotate chimeric.sam with CB:Z: if available
print(f"Annotating {chimeric_sam}...")
with open(chimeric_sam) as fin, open(output_sam, "w") as fout:
    for line in fin:
        if line.startswith("@"):  # header
            fout.write(line)
            continue
        parts = line.strip().split("\t")
        read_name = parts[0]
        if read_name in cb_dict:
            parts.append(f"CB:Z:{cb_dict[read_name]}")
        fout.write("\t".join(parts) + "\n")

print("Done. Output written to:", output_sam)
