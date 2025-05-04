#!/bin/bash

FilePath=$1
mystart=$2
myend=$3
gtffile=$4
mappabilityfile=$5
exonfile=$6
codedir=$7
mkdir -p ${FilePath}/STARMapping
mkdir -p ${FilePath}/ChimericOut
mkdir -p ${FilePath}/Expr/


for ((i=${mystart};i<=${myend};i++))
do
	file=`ls ${FilePath}/STARMapping/${i}/humanChimeric_annotated.out.sam`
	if [[ -n ${file} ]]; then
		python ${codedir}/RmLowMappibility_ChimericRead.py ${file} ${FilePath}/ChimericOut/${i}.sam ${mappabilityfile} 1
	fi
done

bash ${codedir}/Annotate.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${gtffile} ${codedir}
bash ${codedir}/FindFusionSupport.sh ${FilePath}/ChimericOut/ ${mystart} ${myend} ${codedir}
bash ${codedir}/CalcRPKM.sh ${exonfile} ${FilePath}/STARMapping ${FilePath}/Expr/ ${mystart} ${myend}
