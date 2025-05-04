#!/bin/bash

# === CombinePipeline_Predict.sh (updated wrapper including join and Seurat annotate) ===

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
mkdir -p ${FilePath}/FinalResult/

# Step 1: Run MyPredict to generate Prob.txt
python ${codedir}/MyPredict.py ${FilePath}/ChiDist/${prefix}Prob.txt ${weightfile} ${prefix}

# Step 2: Paste ChiDist_middle.txt and Prob.txt to create final ChiDist.txt
paste ${FilePath}/ChiDist/${prefix}ChiDist_middle.txt ${FilePath}/ChiDist/${prefix}Prob.txt > ${FilePath}/ChiDist/${prefix}ChiDist.txt

# Combine all geneanno_cb.sam files
: > "${FilePath}/ChimericOut/geneanno_cb.sam"
for ((i=mystart; i<=myend; i++)); do
  paste "${FilePath}/ChimericOut/${i}_geneanno_cb.sam" >> "${FilePath}/ChimericOut/geneanno_cb.sam"
done

# === Step 3: Join readname_to_genes.txt with each geneanno_cb.sam ===
for ((i=${mystart}; i<=${myend}; i++)); do
  python ${codedir}/JoinGeneannoWithCB.py \
    ${FilePath}/ChiDist/${prefix}readname_to_genes.txt \
    ${FilePath}/ChimericOut/geneanno_cb.sam \
    ${FilePath}/ChiDist/${prefix}readname_to_genes_cb.txt
done

# === Step 4: Run SeuratAnnotate ===
# Use STARsolo outputs per sample in filtered output path

for ((i=${mystart}; i<=${myend}; i++)); do
   sample_path="${FilePath}/STARMapping/${i}/humansolo_Solo.out"  # No need to append /Gene/filtered
   python ${codedir}/seurat_annotate.py $sample_path ${FilePath}/FinalResult/cell_info_${i}.csv
done

# Combine all cell info files (assuming they share headers)
head -n 1 ${FilePath}/FinalResult/cell_info_${mystart}.csv > ${FilePath}/FinalResult/final_cell_info.csv
for ((i=${mystart}; i<=${myend}; i++)); do
  tail -n +2 ${FilePath}/FinalResult/cell_info_${i}.csv >> ${FilePath}/FinalResult/final_cell_info.csv
done

# Step 5: Merge readname_to_genes with cell types
python ${codedir}/MergeGenesWithCellType.py \
  ${FilePath}/ChiDist/readname_to_genes_cb.txt \
  ${FilePath}/FinalResult/final_cell_info.csv \
  ${FilePath}/ChiDist/readname_to_genes_ct.txt

# Step 5: Merge ChiDist with reads
python ${codedir}/MergeChiDistWithRead.py \
  ${FilePath}/ChiDist/${prefix}ChiDist.txt \
  ${FilePath}/ChiDist/readname_to_genes_ct.txt \
  ${FilePath}/ChiDist/${prefix}ChiDist.txt

# Step 6: Filter ChiDist
python ${codedir}/FilterChiDist.py \
  ${FilePath}/ChiDist/${prefix}ChiDist.txt > \
  ${FilePath}/ChiDist/${prefix}ChiDist_filtered.txt