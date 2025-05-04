# Tint Updated for STAR Mapping

#!/bin/bash
filedir=$1
mystart=$2
myend=$3
outdir=$4
genomedir=$5
ncore=$6
whitelist=$7

for ((i=${mystart};i<=${myend};i++)); do
  if [ -f ${filedir}/${i}_1.fastq ]; then
    mkdir -p ${outdir}/${i}

    # === 1. STARsolo ===
    STAR --runThreadN ${ncore} \
      --genomeDir ${genomedir} \
      --readFilesIn ${filedir}/${i}_1.fastq ${filedir}/${i}_2.fastq \
      --soloType CB_UMI_Simple \
      --soloCBstart 1 --soloCBlen 16 \
      --soloUMIstart 17 --soloUMIlen 10 \
      --soloBarcodeReadLength 0 \
      --outSAMattributes NH HI AS nM CB UB \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --soloCBwhitelist ${whitelist} \
      --outFileNamePrefix ${outdir}/${i}/humansolo_

    # === 2. STAR for Chimeric Mapping ===
    STAR --runThreadN ${ncore} \
      --genomeDir ${genomedir} \
      --readFilesIn ${filedir}/${i}_1.fastq ${filedir}/${i}_2.fastq \
      --outSAMtype BAM SortedByCoordinate \
      --chimOutType SeparateSAMold \
      --outSAMunmapped Within KeepPairs \
      --quantMode GeneCounts \
      --outFileNamePrefix ${outdir}/${i}/human \
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
    samtools index ${outdir}/${i}/humansolo_Aligned.sortedByCoord.out.bam

    # === Attach CB tag to Chimeric.out.sam ===
    bam_file="${outdir}/${i}/humansolo_Aligned.sortedByCoord.out.bam"
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