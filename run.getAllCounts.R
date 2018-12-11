args=commandArgs(T)
binCount=185130
sample.list=args[1]
outdir=args[2]
samples=read.table(sample.list,stringsAsFactors=F)[,1]
library(ExomeDepth)
my.counts.matrix <- matrix(ncol=length(samples),nrow=binCount)
colnames(my.counts.matrix) <- samples
for(i in 1:length(samples)){
    my.counts.matrix[,i] <- readRDS(paste0(outdir,"/",samples[i],".my.count.rds"))
}
saveRDS(my.counts.matrix,file=paste0(outdir,"/",'all',".my.counts.rds"))
