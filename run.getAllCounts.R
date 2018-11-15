args=commandArgs(T)
binCount=185130
sample.list=args[1]
samples=read.table(sample.list,stringsAsFactors=F)[,1]
library(ExomeDepth)
my.counts.matrix <- matrix(ncol=length(samples),nrow=binCount)
colnames(my.counts.matrix) <- samples
for(i in 1:length(samples)){
    my.counts.matrix[,i] <- readRDS(paste0(samples[i],".my.count.rds"))
}
saveRDS(my.counts.matrix,file=paste0('all',".my.counts.rds"))
