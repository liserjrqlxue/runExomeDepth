#!/bin/env Rscript
args=commandArgs(T)
sampleID=args[1]
bam=args[2]
gender=args[3]
outdir=args[4]

tag=paste0(".",gender)
suffix=".my.count.rds"

library(ExomeDepth)
source("getBamCount.R")
target <- readRDS("exons.hg19.X.target.rds")
my.count <- getBamCount(
	target = target,
	bam = bam
	)
saveRDS(my.count,file=paste0(outdir,"/",sampleID,tag,suffix))
