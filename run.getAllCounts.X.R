args=commandArgs(T)
binCount=6656
sample.list=args[1]
gender=args[2]
outdir=args[3]

tag=paste0(".",gender)
suffix=".my.count.rds"

samplesInfo<-read.table(sample.list,stringsAsFactors=F)
samples<-samplesInfo[samplesInfo$V3==gender,1]

library(ExomeDepth)
my.counts.matrix <- matrix(ncol=length(samples),nrow=binCount)
colnames(my.counts.matrix) <- samples
for(i in 1:length(samples)){
    message("load ",paste0(outdir,"/",samples[i],tag,suffix))
    my.counts.matrix[,i] <- readRDS(paste0(outdir,"/",samples[i],tag,suffix))
}
saveRDS(my.counts.matrix,file=paste0(outdir,"/",'all',tag,suffix))
