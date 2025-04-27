# Tint Updated for STAR Mapping

#!/bin/bash
# === Arguments ===
filedir=$1          # Input FASTQ directory
mystart=$2          # Start sample index
myend=$3            # End sample index
outdir=$4           # Output directory
genomedir=$5        # STAR genome index directory
ncore=$6            # Number of threads
gtf=$7              # GTF annotation file
whitelist=$8        # Whitelist file
codedir=$9          # Path to the directory containing AddCellBarcodeToChimeric.py

# === Initialize Combined Counts ===
output_counts="${outdir}/combined_counts.txt"
echo -e "GeneID" > $output_counts

# === Main Loop ===
for ((i=${mystart}; i<=${myend}; i++))
do
    if [ -f ${filedir}/${i}_1.fastq ]; then
        mkdir -p ${outdir}/${i}
        echo "Processing sample ${i}..."

        # === STAR Mapping ===
        STAR --runThreadN ${ncore} \
            --genomeDir ${genomedir} \
            --sjdbGTFfile ${gtf} \
            --readFilesIn ${filedir}/${i}_1.fastq ${filedir}/${i}_2.fastq \
            --soloType CB_UMI_Simple \
            --soloCBstart 1 --soloCBlen 16 \
            --soloUMIstart 17 --soloUMIlen 10 \
            --soloBarcodeReadLength 0 \
            --outSAMattributes NH HI AS nM CB UB \
            --outSAMtype BAM SortedByCoordinate \
            --chimOutType SeparateSAMold \
            --outSAMunmapped Within KeepPairs \
            --quantMode GeneCounts \
            --outFileNamePrefix ${outdir}/${i}/ \
            --soloCBwhitelist ${whitelist} \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --chimSegmentReadGapMax 3 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --chimMultimapScoreRange 3 \
            --chimScoreJunctionNonGTAG -4 \
            --chimNonchimScoreDropMin 10 \
            --peOverlapMMp 0.1 \
            > ${outdir}/${i}/STAR_mapping.log 2>&1   # save STAR run log

        # === Update combined gene counts ===
        cut -f 1,2 ${outdir}/${i}/ReadsPerGene.out.tab | sed '1d' >> $output_counts

        # === Attach CB tag to Chimeric.out.sam ===
        bam_file="${outdir}/${i}/Aligned.sortedByCoord.out.bam"
        chimeric_sam="${outdir}/${i}/Chimeric.out.sam"
        output_sam="${outdir}/${i}/Chimeric_annotated.out.sam"

        echo "Annotating chimeric reads for sample ${i}..."
        python ${codedir}/AddCellBarcodeToChimeric.py "$bam_file" "$chimeric_sam" "$output_sam" \
            > ${outdir}/${i}/barcode_annotation.log 2>&1  # save barcode annotation log

        echo "Sample ${i} finished!"
    else
        echo "Warning: ${filedir}/${i}_1.fastq not found, skipping sample ${i}."
    fi
done

echo "All samples processed!"
