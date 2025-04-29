# Tint Updated to support one extra column

from __future__ import print_function
from __future__ import division
import sys
import numpy

# ***** readme *****
# This code filters the fusion candidates that have fewer than 3 supporters.

def MakeString(mylist, sep='\t'):
    string = ''
    for i in range(len(mylist) - 1):
        string += str(mylist[i]) + sep
    string += mylist[-1]
    return string

bigcount = 0
ChiDistFile = open(sys.argv[1])
templines = []
linedic = {}

for line in ChiDistFile.readlines():
    if line[0] == '#':
        continue
    try:
        info = line.strip().split('\t')
        
        # Extract split_scores (info[3])
        cc = info[3]
        count = cc.strip('[]').replace("'", "").split(', ')
        
        # Clean empty or invalid counts
        clean_count = []
        for item in count:
            try:
                if int(item) >= 1:
                    clean_count.append(item)
            except:
                continue
        
        mysum = sum(int(x) for x in clean_count)
        
        if len(clean_count) <= 1:
            continue
        if mysum / len(clean_count) > 5:
            continue
        
        # Update fields
        info[3] = '[' + ', '.join(clean_count) + ']'
        info[2] = str(len(clean_count))
        
        # Reconstruct line
        templines.append(MakeString(info))
        
    except Exception as e:
        pass  # Silently skip bad lines

ChiDistFile.close()

for line in templines:
    info = line.strip().split('\t')
    score = len(info[3])  # Still use split score string length for ranking
    linedic[line] = score

numrecord = len(linedic)

for key in sorted(linedic, key=linedic.__getitem__, reverse=True):
    if bigcount < min(10.0, numpy.floor(numrecord / 50)):
        bigcount += 1
    else:
        print(key, end='')
