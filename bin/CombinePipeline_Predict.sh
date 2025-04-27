# Tint updated for Seurat Annotate

#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
prefix=$4
weightfile=$5
hg19file=$6
gtf=$7
codedir=$8

if [ "${prefix}" = "." ]
then
	prefix=""
fi

mkdir -p ${FilePath}/ChiDist/

# Step 1: Run MyPredict to generate Prob.txt
python ${codedir}/MyPredict.py ${FilePath}/ChiDist/${prefix}Prob.txt ${weightfile} ${prefix}

# Step 2: Paste ChiDist_middle.txt and Prob.txt to create final ChiDist.txt
paste ${FilePath}/ChiDist/${prefix}ChiDist_middle.txt ${FilePath}/ChiDist/${prefix}Prob.txt > ${FilePath}/ChiDist/${prefix}ChiDist.txt

# Step 3: Run SeuratAnnotate (NEW)
# Assuming your STARsolo filtered outputs are stored at ${FilePath}/STARsoloOutput/
Rscript ${codedir}/SeuratAnnotate.R ${FilePath}/STARsoloOutput/ ${FilePath}/FinalResult/final_cell_info.csv
python ${codedir}/MergeChiDistWithCellType.py ${FilePath}/ChiDist/${prefix}ChiDist.txt ${FilePath}/FinalResult/final_cell_info.csv ${FilePath}/ChiDist/${prefix}ChiDist.txt

# Step 4: Filter ChiDist
python ${codedir}/FilterChiDist.py ${FilePath}/ChiDist/${prefix}ChiDist.txt > ${FilePath}/ChiDist/${prefix}ChiDist_filtered.txt
