#!/bin/env Rscript
args=commandArgs(T)
ref=args[1]
sampleID=args[2]
bam=args[3]
tag=args[4]
outdir=args[5]

suffix=".my.count.rds"
suffix=paste0(".",tag,suffix)

library(ExomeDepth)
if(length(args)>5){
	sourceDir=args[6]
	source(paste(sourceDir,"getBamCount.R",sep="/"))
	if(tag=="A"){
		target <- readRDS(paste0(sourceDir,"/exons.",ref,".A.target.rds"))
	}else{
		target <- readRDS(paste0(sourceDir,"/exons.",ref,".X.target.rds"))
	}
}else{
	source("getBamCount.R")
	if(tag=="A"){
		target <- readRDS(paste0("exons.",ref,".A.target.rds"))
	}else{
		target <- readRDS(paste0("exons.",ref,".X.target.rds"))
	}
}
my.count <- getBamCount(
	target = target,
	bam = bam
	)
saveRDS(my.count,file=paste0(outdir,"/",sampleID,suffix))
