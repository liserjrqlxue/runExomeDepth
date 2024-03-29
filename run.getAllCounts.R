args=commandArgs(T)
sample.list=args[1]
tag=args[2]
outdir=args[3]

suffix=".my.count.rds"
suffix=paste0(".",tag,suffix)

if(tag=="A"){
	binCount=185130
	samples<-read.table(sample.list,stringsAsFactors=F)[,1]
}else{
	binCount=6656
	samplesInfo<-read.table(sample.list,stringsAsFactors=F,colClasses = c("character"))
	samples<-samplesInfo[samplesInfo$V2==tag,1]
}

library(ExomeDepth)
my.counts.matrix <- matrix(ncol=length(samples),nrow=binCount)
colnames(my.counts.matrix) <- samples
if(length(samples)>0){
	for(i in 1:length(samples)){
	    message("load ",paste0(outdir,"/",samples[i],suffix))
	    my.counts.matrix[,i] <- readRDS(paste0(outdir,"/",samples[i],suffix))
	}
}else{
    message("no ",paste0(suffix)," to load!")
}

if(length(args)>3){
	controlDir=args[4]
	controlMatrix<-readRDS(paste0(controlDir,"/all",suffix))
	my.counts.matrix<-cbind(my.counts.matrix,controlMatrix)
}

saveRDS(my.counts.matrix,file=paste0(outdir,"/",'all',suffix))
