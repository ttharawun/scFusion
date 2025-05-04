#!/usr/bin/env python3

# Tint Updated: Merge ChiDist.txt with Seurat CellType using CB (column 21)

import sys
import pandas as pd

# Arguments
chifile = sys.argv[1]       # Path to ChiDist.txt
cellinfofile = sys.argv[2]  # Path to final_cell_info.csv
outputfile = sys.argv[3]    # Path to write updated ChiDist_with_celltype.txt

# Step 1: Read ChiDist file (tab-delimited, no header)
chi = pd.read_csv(chifile, sep='\t', header=None, dtype=str)

# Validate expected structure (CB should be at index 21)
if chi.shape[1] < 22:
    sys.exit("Error: ChiDist file must have at least 22 columns (including CB at col 21).")

# Step 2: Read Seurat cell type info
cellinfo = pd.read_csv(cellinfofile)
if 'Barcode' not in cellinfo.columns or 'CellType' not in cellinfo.columns:
    sys.exit("Error: final_cell_info.csv must have columns 'Barcode' and 'CellType'.")

# Step 3: Create mapping dict and annotate
celltype_dict = dict(zip(cellinfo['Barcode'], cellinfo['CellType']))

def map_celltype(barcode):
    return celltype_dict.get(barcode, "Unknown")

# Step 4: Append CellType column after CB (col 21)
chi[22] = chi[21].apply(map_celltype)

# Step 5: Save to output
chi.to_csv(outputfile, sep='\t', header=False, index=False)

print(f"✅ Merged {len(chi)} rows. Output saved to: {outputfile}")
