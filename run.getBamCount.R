#!/bin/env Rscript
args=commandArgs(T)
sampleID=args[1]
bam=args[2]
library(ExomeDepth)
source("getBamCount.R")
target <- readRDS("exons.hg19.target.rds")
my.count <- getBamCount(
	target = target,
	bam = bam
	)
saveRDS(my.count,file=paste0(sampleID,".my.count.rds"))
