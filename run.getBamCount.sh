#!/bin/env bash
set -euo pipefail
Bin=$(dirname $(readlink -f "$0"))
sampleID=$1
bam=$2
gender=$3
outDir=$4

# input  : $bam
# output : $outDir/$sampleID.A.my.count.rds $outDir/$sampleID.$gender.my.count.rds

mkdir -p $outDir

Rscript $Bin/run.getBamCount.R $sampleID $bam $gender $outDir $Bin
