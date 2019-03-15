# ExomeDepth for WES CNV calling

## Base on https://github.com/vplagnol/ExomeDepth.git

## simple run
```bash
Rscript run.R bam.list bai.list
```

## speed up
1. create Bins
  * Input:  
`ref` as hg19 reference file
  * Output:  
`exons.hg19.A.target.rds`  
`exons.hg19.A.rdata.rds`
`exons.hg19.X.target.rds`  
`exons.hg19.X.rdata.rds`
  * Script:  
```R
library(ExomeDepth)
data(exons.hg19)
data(exons.hg19.X)
source("getBins.R")
ref='hg19.fa'
getBins(
    bed.frame = exons.hg19,
    include.chr = T,
    referenceFasta = ref,
    prefix="exons.hg19.A"
    )
getBins(
    bed.frame = exons.hg19.X,
    include.chr = T,
    referenceFasta = ref,
    prefix="exons.hg19.X"
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
source("getBamCount.R")
target <- readRDS("exons.hg19.target.rds")
my.count <- getBamCount(
    target = target,
    bam = bam
    )
saveRDS(my.count,file=paste0(sampleID,".my.count.rds"))
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
#sample.list="sample.list"
samples=read.table(sample.list,stringsAsFactors=F)[,1]
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
samples <- colnames(my.counts.matrix)
for(i in 1:length(samples)){
    if(sampleID!=samples[i]){
        next
    }
    my.choice <- select.reference.set(
        test.counts = my.counts.matrix[,i],
        reference.counts = my.counts.matrix[,-i],
        bin.length = (exons.hg19$end - exons.hg19$start)/1000,
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
        chromosome = exons.hg19$chromosome,
        start = exons.hg19$start+1,
        end = exons.hg19$end,
        name = exons.hg19$name
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
        reference.annotation = exons.hg19.GRanges,
        min.overlap = 0.0001,
        column.name = 'Conrad.hg19'
        )
    output.file <- paste0(sampleID,'.CNV.calls.tsv')
    write.table(file=output.file,x=all.exons@CNV.calls,row.names=F,sep="\t",quote=F)
    saveRDS(all.exons,file=paste0(sampleID,'.all.exons.rds'))
}
```

## example
1. `Rscript run.getBamCount.R SampleID bam.file` for each sample  
input:`SampleID` as prefix and `bam.file` as bam path  
output:`SampleID.my.count.rds`  
time:15mins for each  
2. `Rscript run.getAllCounts.R sample.list`  
input:`sample.list` as batch list of samples  
output:`all.my.counts.rds`  
time:12s  
3. `Rscript run.CallCNVs.R`  
input:`all.my.counts.rds`  
output:`SampleID.CNV.calls.tsv` and `SampleID.all.exons.rds` for each sample  
time:12mins for total 12 test samples  
4. `Rscript run.getCNVs.R SampleID`  
input:`all.my.counts.rds` and selected sample `SampleID`  
output:`SampleID.CNV.calls.tsv` and `SampleID.all.exons.rds` for selected sample  
time:1mins
