# Tint added for Artifact Scoring step

# JoinGeneannoWithCB.py
import sys

gene_read_file = sys.argv[1]   # readname_to_genes.txt
sam_file        = sys.argv[2]   # geneanno_cb.sam
output_file     = sys.argv[3]   # e.g. readname_to_genes_CB.txt

# Step 1: Build a dict mapping (gene, chr, pos) → CB
cb_coord = {}
with open(sam_file) as sam:
    for line in sam:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        gene = fields[0]
        rname = fields[3]
        pos   = fields[4]
        # find the CB:Z: tag
        cb = "-"
        for tag in fields[11:]:
            if tag.startswith("CB:Z:"):
                cb = tag[5:]
                break
        key = (gene, rname, pos)
        # prefer a real barcode over "-"
        if key not in cb_coord or (cb != "-" and cb):
            cb_coord[key] = cb

# Step 2: Annotate each fusion read by checking geneA site first, then geneB
with open(gene_read_file) as infile, open(output_file, "w") as out:
    header = infile.readline().rstrip("\n")
    out.write(header + "\tCB\n")

    for line in infile:
        parts = line.rstrip("\n").split("\t")
        # if missing coords, just write line + “\t–” and continue
        if len(parts) < 7:
            out.write(line.rstrip("\n") + "\t-\n")
            continue
        read_name, geneA, geneB, chr1, chr2, pos1, pos2 = parts

        # lookup CB: try geneA site first, then geneB
        cb = cb_coord.get((geneA, chr1, pos1), "-")
        if cb == "-" or not cb:
            cb = cb_coord.get((geneB, chr2, pos2), "-")

        out.write(line.rstrip("\n") + f"\t{cb}\n")
