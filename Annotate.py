from __future__ import print_function
import sys
from intervaltree import Interval, IntervalTree

# --- Parse GTF file into interval trees ---
def parse_gtf(gtf_file):
    gene_trees = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene':
                chrom = fields[0]
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom  # Ensure 'chr' prefix

                start = int(fields[3])
                end = int(fields[4])
                info = fields[8]

                gene_name = None
                for entry in info.split(';'):
                    if 'gene_name' in entry:
                        gene_name = entry.split('"')[1]
                        break

                if gene_name:
                    if chrom not in gene_trees:
                        gene_trees[chrom] = IntervalTree()
                    gene_trees[chrom][start:end+1] = gene_name  # end+1 because intervaltree is [start, end)
    return gene_trees

# --- Find gene by chromosome and position ---
def find_gene(gene_trees, chrom, pos):
    if chrom in gene_trees:
        hits = gene_trees[chrom][pos]
        if hits:
            # If multiple genes overlap, just pick one (can improve later)
            return list(hits)[0].data
    return None

# --- Main Program ---
if __name__ == "__main__":
    input_sam = sys.argv[1]
    gtf_file = sys.argv[2]
    output_sam = input_sam[:-4] + '_geneanno.sam'

    # Step 1: Parse GTF into interval trees
    print("Parsing GTF...")
    gene_trees = parse_gtf(gtf_file)
    print("Done parsing GTF.")

    # Step 2: Process SAM file
    with open(input_sam) as fin, open(output_sam, 'w') as fout:
        for line in fin:
            if line.startswith('@'):
                fout.write(line)
                continue
            info = line.strip().split('\t')
            chrom = info[2]
            pos = int(info[3])

            genename = find_gene(gene_trees, chrom, pos)
            if not genename:
                genename = 'NA'

            fout.write(genename + '\t' + line)
