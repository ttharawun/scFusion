# Tint Updated for STAR Mapping

#!/bin/bash
# === Arguments ===
filedir=$1          # Input FASTQ directory
mystart=$2          # Start sample index
myend=$3            # End sample index
outdir=$4           # Output directory
genomedir=$5        # STAR genome index directory
gtf=$6              # GTF annotation file
whitelist=$7        # Whitelist file
ncore=$8            # Number of threads

# === Automatically set code directory ===
codedir=$(dirname "$0")   # Path to bin/ folder

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
            --outFileNamePrefix ${outdir}/${i}/human \
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
            --peOverlapMMp 0.1

	# === Index BAM ===
        samtools index ${outdir}/${i}/humanAligned.sortedByCoord.out.bam

        # === Attach CB tag to Chimeric.out.sam ===
        bam_file="${outdir}/${i}/humanAligned.sortedByCoord.out.bam"
        chimeric_sam="${outdir}/${i}/humanChimeric.out.sam"
        output_sam="${outdir}/${i}/humanChimeric_annotated.out.sam"

        echo "Annotating chimeric reads for sample ${i}..."
        python ${codedir}/AddCellBarcodeToChimeric.py "$bam_file" "$chimeric_sam" "$output_sam"

        echo "Sample ${i} finished!"
    else
        echo "Warning: ${filedir}/${i}_1.fastq not found, skipping sample ${i}."
    fi
done

echo "All samples processed!"

