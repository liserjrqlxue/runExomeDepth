#!/bin/bash
set -euo pipefail
Bin=$(dirname $(readlink -f "$0"))
sampleID=$1
gender=$2
outDir=$3

# input  : $outDir/all.$gender.my.count.rds
# output : $outDir/$sampleID.$gender.CNV.calls.tsv $outDir/$sampleID.$gender.all.exons.rds

Rscript $Bin/run.getCNVs.R $sampleID $gender $outDir
