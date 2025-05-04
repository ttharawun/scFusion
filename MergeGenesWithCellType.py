# Tint added for merging genes with cell types 

#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: python MergeGenesWithCellType.py "
            "<readname_to_genes.txt> <cell_info.csv> <output.txt>\n"
        )
        sys.exit(1)

    genes_path   = sys.argv[1]
    cells_path   = sys.argv[2]
    out_path     = sys.argv[3]

    # load tables
    genes = pd.read_csv(genes_path, sep="\t", dtype=str)
    cells = pd.read_csv(cells_path, dtype=str)

    # sanity checks
    if "CB" not in genes.columns:
        sys.stderr.write("ERROR: 'CB' column not found in readname_to_genes.txt\n")
        sys.exit(1)
    if not {"Barcode", "CellType"}.issubset(cells.columns):
        sys.stderr.write("ERROR: 'Barcode' and/or 'CellType' missing in cell_info.csv\n")
        sys.exit(1)

    # merge on CB ↔ Barcode
    merged = genes.merge(
        cells[["Barcode","CellType"]],
        how="left",
        left_on="CB",
        right_on="Barcode"
    )

    # drop the extra Barcode, rename CellType→CT, fill missing
    merged.drop(columns="Barcode", inplace=True)
    merged.rename(columns={"CellType": "CT"}, inplace=True)
    merged["CT"] = merged["CT"].fillna("-")

    # write out tab-delimited
    merged.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
