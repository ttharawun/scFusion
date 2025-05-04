from __future__ import print_function
import sys
from intervaltree import IntervalTree

# --- Parse GTF file into interval trees ---
def parse_gtf(gtf_file):
    gene_trees = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] != 'gene':
                continue

            chrom = fields[0]
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom  # Enforce 'chr' prefix

            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            attrs = fields[8]

            gene_name = ""
            for attr in attrs.strip().split(';'):
                if 'gene_name' in attr:
                    gene_name = attr.strip().split('"')[1]
                    break

            if gene_name:
                if chrom not in gene_trees:
                    gene_trees[chrom] = IntervalTree()
                gene_trees[chrom][start:end] = gene_name
    return gene_trees

# --- Find gene by chromosome and position with ±20bp fallback ---
def find_gene(gene_trees, chrom, pos):
    hits = set()
    for delta in [0, 20, -20]:
        query_pos = pos + delta
        if chrom in gene_trees:
            found = gene_trees[chrom][query_pos]
            hits.update(g.data for g in found)
        if hits:
            break
    return ";".join(hits) if hits else "NA"

# --- Main Execution ---
if __name__ == "__main__":
    input_sam = sys.argv[1]
    gtf_file = sys.argv[2]
    output_nocb = input_sam[:-4] + "_geneanno.sam"
    output_withcb = input_sam[:-4] + "_geneanno_cb.sam"

    print("Parsing GTF file...")
    gene_trees = parse_gtf(gtf_file)
    print("Finished parsing GTF.")

    with open(input_sam) as fin, \
         open(output_nocb, 'w') as fout_nocb, \
         open(output_withcb, 'w') as fout_withcb:

        for line in fin:
            if line.startswith('@'):
                fout_nocb.write(line)
                fout_withcb.write(line)
                continue

            fields = line.strip().split('\t')
            chrom = fields[2]
            pos = int(fields[3]) - 1  # SAM is 1-based

            genename = find_gene(gene_trees, chrom, pos)

            # For _geneanno_cb.sam (keep original CB:Z: tag)
            fout_withcb.write(genename + '\t' + line)

            # For _geneanno.sam (remove CB:Z:* field)
            core = fields[:11]
            optional = [f for f in fields[11:] if not f.startswith("CB:Z:")]
            cleaned_line = '\t'.join(core + optional)
            fout_nocb.write(genename + '\t' + cleaned_line + '\n')
