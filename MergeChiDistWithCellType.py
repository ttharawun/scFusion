# Tint Added to merge cell type and fusion report

#!/usr/bin/env python

import sys
import pandas as pd

# Arguments
chifile = sys.argv[1]     # Path to ChiDist.txt
cellinfofile = sys.argv[2] # Path to final_cell_info.csv
outputfile = sys.argv[3]   # Path to write updated ChiDist_with_celltype.txt

# Read ChiDist
chi = pd.read_csv(chifile, sep='\t', header=None, dtype=str)

# Read Cell Info
cellinfo = pd.read_csv(cellinfofile)

# Build barcode → cell type dictionary
celltype_dict = dict(zip(cellinfo['Barcode'], cellinfo['CellType']))

# Map CB (column 18) to cell type
def map_celltype(barcode):
    return celltype_dict.get(barcode, "Unknown")  # If not found, label as "Unknown"

# Add new column
chi['CellType'] = chi[18].apply(map_celltype)

# Save merged ChiDist
chi.to_csv(outputfile, sep='\t', header=False, index=False)
