#!/bin/env Rscript
args=commandArgs(T)
sampleID=args[1]
bam=args[2]
outdir=args[3]
library(ExomeDepth)
if(length(args)>3){
	sourceDir=args[4]
	source(paste(sourceDir,"getBamCount.R",sep="/"))
	target <- readRDS(paste(sourceDir,"exons.hg19.target.rds",sep="/"))
}else{
	source("getBamCount.R")
	target <- readRDS("exons.hg19.target.rds")
}
my.count <- getBamCount(
	target = target,
	bam = bam
	)
saveRDS(my.count,file=paste0(outdir,"/",sampleID,".my.count.rds"))
