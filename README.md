# ExomeDepth for WES CNV calling

## Base on https://github.com/vplagnol/ExomeDepth.git

## simple run
```bash
Rscript run.R bam.list bai.list
```

## TO-DO:
### speed up
1. create Bins
  * Input:  
`ref` as hg19 reference file
  * Output:  
`exons.hg19.target.rds`  
`exons.hg19.rdata.rds`
  * Script:  
```R
library(ExomeDepth)
data(exons.hg19)
source(getBins.R)
ref='hg19.fa'
getBins(
    bed.frame = exons.hg19,
    include.chr = T,
    referenceFasta = ref
    )
```

2. run getBamCounts parallel for each sample
  * Input:  
`sampleID` as prefix  
`bam` as bam file path
  * Output:  
`SampleID.my.count.rds`
  * Script:  
```R
args=commandArgs(T)
sampleID=args[1]
bam=args[2]
library(ExomeDepth)
data(exons.hg19)
source(getBamCount.R)
target <- readRDS(exons.hg19.target.rds)
my.count <- getBamCount(
    target = target,
    bam = bam.file
    )
saveRDS(my.count,file=paste0(sampleID,".my.count.rds))
```

3. create main matrix of read count data of all sample
  * Input:  
`sample.list` as `sampleID` list
  * Output:  
`all.my.counts.rds` as all count data in a matrix
  * Script:  
```R
args=commandArgs(T)
sample.list=args[1]
samples=read.table(sample.list,stringAsFactors=F)[,1]
library(ExomeDepth)
my.bins <- readRDS('exons.hg19.rdata.rds')
my.counts.matrix <- matrix(ncol=length(samples),nrow=nrow(my.bins))
colnames(my.counts.matrix) <- samples
for(i in 1:length(samples)){
    my.counts.matrix[,i] <- readRDS(paste0(samples[i],".my.count.rds"))
}
saveRDS(my.counts.matrix,file=paste0('all',".my.counts.rds"))
```

4. run CallCNVs for each sample
  * Input:  
`sampleID` as selected sample to call CNVs
  * Output:  
`sampleID.CNV.calls.tsv` as CNV calls result
`sampleID.all.exons.rds` as rds data of CNV calls
  * Script:  
```R
args=commandArgs(T)
sampleID=args[1]
library(ExomeDepth)
data(exons.hg19)
data(Conrad.hg19)
exons.hg19.GRanges <- GenomicRanges::GRanges(
    seqnames = exons.hg19$chromosome,
    IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
    names = exons.hg19$name
    )
my.counts.matrix <- readRDS(file=paste0('all',".my.counts.rds"))
my.bins <- readRDS('exons.hg19.rdata.rds')
my.bins.dafr <- as(my.bins[,colnames(my.bins)],'data.frame')
my.bins.dafr$chromosome <-
    gsub(
        as.character(my.bins.dafr$space),
        pattern = 'chr',
        replacement = ''
        )
samples <- colnames(my.counts.matrix)
for(i in 1:length(samples)){
    if(sampleID!=samples[i]){
        next
    }
    my.choice <- select.reference.set(
        test.counts = my.counts.matrix[,i],
        reference.counts = my.counts.matrix[,-i],
        bin.length = (my.bins.dafr$end - my.bins.dafr$start)/1000,
        n.bins.reduced = 10000
        )
    my.reference.selected <- apply(
        X = my.counts.matrix[,my.choice$reference.choice,drop=FALSE],
        MAR = 1,
        FUN = sum
        )
    message('Now creating the ExomeDepth object')
    all.exons <- new(
        'ExomeDepth',
        test = my.counts.matrix[,i],
        reference = my.reference.selected,
        formula = 'cbind(test, reference) ~ 1'
        )
######## Now call the CNVs
    all.exons <- CallCNVs(
        x = all.exons,
        transition.probability = 10^-4,
        chromosome = my.counts.dafr$chromosome,
        start = my.counts.dafr$start,
        end = my.counts.dafr$end,
        name = my.counts.dafr$names
        )
######## Now annotate the ExomeDepth object
    all.exons <- AnnotateExtra(
        x = all.exons,
        reference.annotation = Conrad.hg19.common.CNVs,
        min.overlap = 0.5,
        column.name = 'Conrad.hg19'
        )
    all.exons <- AnnotateExtra(
        x = all.exons,
        reference.annotation = exons.gh19,GRanges,
        min.overlap = 0.0001,
        column.name = 'Conrad.hg19'
        )
    output.file <- paste0(sampleID,'.CNV.calls.tsv')
    write.table(file=output.file,x=all.exons@CNV.calls,row.names=F,sep="\t",quote=F)
    saveRDS(all.exons,file=paste0(sampleID,'.all.exons.rds'))
}
```
