#!/bin/env Rscript
args=commandArgs(T)
sampleID=args[1]
bam=args[2]
tag=args[3]
outdir=args[4]

suffix=".my.count.rds"
suffix=paste0(".",tag,suffix)

library(ExomeDepth)
if(length(args)>4){
	sourceDir=args[5]
	source(paste(sourceDir,"getBamCount.R",sep="/"))
	if(tag=="A"){
		target <- readRDS(paste(sourceDir,"exons.hg19.A.target.rds",sep="/"))
	}else{
		target <- readRDS(paste(sourceDir,"exons.hg19.X.target.rds",sep="/"))
	}
}else{
	source("getBamCount.R")
	if(tag=="A"){
		target <- readRDS("exons.hg19.A.target.rds")
	}else{
		target <- readRDS("exons.hg19.X.target.rds")
	}
}
my.count <- getBamCount(
	target = target,
	bam = bam
	)
saveRDS(my.count,file=paste0(outdir,"/",sampleID,suffix))
