# Tint Updated to include cell barcodes

#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys
import math
import random

# ***** readme *****
# This code extracts chimeric read from sam file for training, with pos and direction
# The input is *ChiDist_middle.txt

random.seed(1122)
infile = open(sys.argv[1])
lines = infile.readlines()
totallines = len(lines)
uselines = random.sample(
    lines[int(totallines/10):],
    math.floor(min(int(totallines * 0.4), max(15000, totallines/9)))
)

for line in uselines:
    info = line.rstrip().split('\t')
    gene1 = info[0]
    gene2 = info[1]
    if gene1.startswith('IG') or gene2.startswith('IG') \
       or gene1.startswith('TRA') or gene2.startswith('TRA'):
        continue

    # extract training fields
    read     = info[-5]
    splitpos = info[-4]
    dir1     = info[-3]  # original “direct1”
    dir2     = info[-2]  # original “direct2”
    cb       = info[-1]  # your cell barcode

    # print four‐field training line + CB as fifth column
    print(
        read     + '\t' +
        splitpos + '\t' +
        f"{info[6]}:{info[8]}:{dir1}" + '\t' +
        f"{info[7]}:{info[9]}:{dir2}" + '\t' +
        cb
    )

infile.close()
