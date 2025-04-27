# Tint Updated to include cell barcodes

from __future__ import print_function
import sys

# Extract specific tag from optional SAM fields
def extract_tag(info, tag_prefix):
    for field in info:
        if field.startswith(tag_prefix):
            return field[len(tag_prefix):]
    return "NA"

def ReverseComplement(str):
    return str[::-1].replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

def SolveClip(str, readlen):
    num = []
    alp = []
    res = []
    start = 0
    scount = 0
    Ndel = 0
    length = len(str)
    lastalp = -1
    for i in range(length):
        if str[i].isalpha():
            alp.append(str[i])
            num.append(int(str[lastalp + 1:i]))
            lastalp = i
    i = 0
    while i < len(alp):
        if alp[i] == 'N':
            for j in range(i):
                if alp[j] == 'S':
                    scount = 1
                    break
            del alp[i]
            Ndel += num[i]
            del num[i]
            continue
        i += 1
    lastm = -1
    i = 0
    while i < len(alp):
        if alp[i] == 'M':
            if lastm == -1:
                lastm = i
                continue
            if lastm == i - 1:
                num[i - 1] += num[i]
                del alp[i]
                del num[i]
                continue
            lastm = i
        i += 1
    numsum = [num[0]]
    for i in range(1, len(num)):
        numsum.append(numsum[-1] + num[i])
    for i in range(len(alp)):
        if alp[i] == 'M':
            if i == 0:
                res = [numsum[0]]
            elif i == len(alp) - 1:
                res.append(numsum[i - 1])
            else:
                res.append(numsum[i - 1])
                res.append(numsum[i])
    for i in range(len(alp)):
        if alp[i] == 'M':
            if i > 0:
                start = numsum[i - 1]
            break
    return [res, alp, start, 0, Ndel, scount]

def TakeoutFusionSupport(lines):
    if len(lines) == 3:
        info1 = lines[0].split('\t')
        info2 = lines[1].split('\t')
        info3 = lines[2].split('\t')
        gene1 = info1[0]
        gene2 = info2[0]
        gene3 = info3[0]
        if (gene1 != gene2 or gene2 != gene3) and (gene1 == gene2 or gene2 == gene3 or gene1 == gene3) and (gene1 != '' and gene2 != '' and gene3 != ''):
            if info1[10] == info2[10] or info1[10] == ReverseComplement(info2[10]):
                readinfo1 = info1
                readinfo2 = info2
            elif info1[10] == info3[10] or info1[10] == ReverseComplement(info3[10]):
                readinfo1 = info1
                readinfo2 = info3
            elif info2[10] == info3[10] or info2[10] == ReverseComplement(info3[10]):
                readinfo1 = info3
                readinfo2 = info2
            else:
                print('Three reads are different!')
                for ll in lines:
                    print(ll, end='')
                return thisname
            if readinfo1[0] != readinfo2[0]:
                clip1 = readinfo1[6]
                clip2 = readinfo2[6]
                readlen = len(readinfo1[10])
                clipsplit1 = SolveClip(clip1, readlen)
                clipsplit2 = SolveClip(clip2, readlen)
                cc = False
                splitpnt1 = -1
                splitpnt2 = -1
                for i in range(len(clipsplit1[0])):
                    for j in range(len(clipsplit2[0])):
                        if clipsplit2[0][j] == clipsplit1[0][i] or clipsplit2[0][j] + clipsplit1[0][i] == readlen:
                            splitpnt1 = clipsplit1[0][i]
                            splitpnt2 = clipsplit2[0][j]
                            cc = True
                            break
                    if cc:
                        break
                if readlen - 0.5 > splitpnt1 > -0.5 and readlen - 0.5 > splitpnt2 > 0.5:
                    if clipsplit1[1][i] == 'S':
                        brkpnt1 = int(readinfo1[4])
                        direct1 = 1
                    else:
                        brkpnt1 = splitpnt1 - clipsplit1[2] + int(readinfo1[4]) + clipsplit1[4] - 1
                        direct1 = -1
                    if clipsplit2[1][j] == 'S':
                        brkpnt2 = int(readinfo2[4])
                        direct2 = 1
                    else:
                        brkpnt2 = splitpnt2 - clipsplit2[2] + int(readinfo2[4]) + clipsplit2[4] - 1
                        direct2 = -1
                    cgene1 = readinfo1[0]
                    cgene2 = readinfo2[0]
                    chromo1 = readinfo1[3]
                    chromo2 = readinfo2[3]
                    if not chromo1.startswith('chr'):
                        chromo1 = 'chr' + chromo1
                    if not chromo2.startswith('chr'):
                        chromo2 = 'chr' + chromo2
                    cb_tag1 = extract_tag(readinfo1, "CB:Z:")
                    cb_tag2 = extract_tag(readinfo2, "CB:Z:")
                    cb_set = set([cb_tag1, cb_tag2])
                    genelist1 = cgene1.split(';')
                    genelist2 = cgene2.split(';')
                    for gene1 in genelist1:
                        for gene2 in genelist2:
                            fusion_key = gene1 + '\t' + gene2
                            reverse_key = gene2 + '\t' + gene1
                            entry = [brkpnt1, brkpnt2, splitpnt1, clipsplit1[2], clipsplit2[2], direct1, direct2]
                            if fusion_key in geneset:
                                geneset[fusion_key][1] += 1
                                geneset[fusion_key][2].append(entry)
                                geneset[fusion_key][4].update(cb_set)
                            elif reverse_key in geneset:
                                geneset[reverse_key][1] += 1
                                geneset[reverse_key][2].append(entry[::-1])
                                geneset[reverse_key][4].update(cb_set)
                            else:
                                geneset[fusion_key] = [0, 1, [entry], [chromo1, chromo2], set(cb_set)]
    return ''

geneset = {}
infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')
lastname = ''
lines = []
badreadname = ''
thisname = 'hgfhfrjfjzjbest1122gh'

with infile:
    for line in infile:
        if len(line) > 0 and line[0] == '#':
            continue
        info = line.split('\t')
        if len(info) > 1:
            thisname = info[1]
        if thisname == badreadname:
            continue
        if thisname == lastname:
            lines.append(line)
        else:
            aa = TakeoutFusionSupport(lines)
            lastname = aa if len(aa) > 1 else thisname
            lines = [line]
        if len(info[0]) <= 1:
            badreadname = info[1]
            lastname = ''
            lines = []
aa = TakeoutFusionSupport(lines)

for key in geneset:
    outline = key + '\t' + str(geneset[key][0]) + '\t' + str(geneset[key][1]) + '\t' + geneset[key][3][0] + '\t' + geneset[key][3][1] + '\t'
    posstr = []
    for i in geneset[key][2]:
        posstr.append(str(i[0]) + ',' + str(i[1]) + '+' + str(i[2]) + '+' + str(i[3]) + ',' + str(i[4]) + ';')
    outline += ''.join(set(posstr))
    outline += '\t' + ';'.join(sorted(geneset[key][4]))  # Add CB tags
    outfile.write(outline + '\n')

infile.close()
outfile.close()
