#!/usr/bin/env python3
import sys
import pandas as pd

def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: python MergeChiDistWithRead.py "
            "<ChiDist.txt> <readname_to_genes_CT.txt> <output.txt>\n"
        )
        sys.exit(1)

    chi_path   = sys.argv[1]
    annot_path = sys.argv[2]
    out_path   = sys.argv[3]

    # 1) Read ChiDist without header
    chi = pd.read_csv(chi_path, sep="\t", header=None, dtype=str)
    # assign names so we can merge on the first 10 cols
    n_cols = chi.shape[1]
    base = ['geneA','geneB','_2','_3','_4','_5','chr1','chr2','pos1','pos2']
    extra = [f'col{i}' for i in range(10, n_cols)]
    chi.columns = base + extra

    # 2) Read annotations (must have CB & CT)
    ann = pd.read_csv(annot_path, sep="\t", dtype=str)
    for col in ['geneA','geneB','chr1','chr2','pos1','pos2','CB','CT']:
        if col not in ann.columns:
            sys.stderr.write(f"ERROR: missing column {col} in {annot_path}\n")
            sys.exit(1)

    # 3) Merge on the six keys
    merged = chi.merge(
        ann[['geneA','geneB','chr1','chr2','pos1','pos2','CB','CT']],
        on=['geneA','geneB','chr1','chr2','pos1','pos2'],
        how='left'
    )

    # 4) Fill missing CB/CT
    merged['CB'] = merged['CB'].fillna('-')
    merged['CT'] = merged['CT'].fillna('-')

    # 5) Write out WITHOUT header or index
    merged.to_csv(out_path, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
