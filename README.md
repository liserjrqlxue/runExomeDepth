# ExomeDepth for WES CNV calling

## simple run
```bash
Rscript run.R bam.list bai.list
```

## TO-DO:
### speed up
1. create Bins
  * Script
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
  * input:  
`ref` as hg19 reference file
  * output:  
`exons.hg19.target.rds`  
`exons.hg19.rdata.rds`
2. run getBamCounts parallel for each sample  
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
  * input:
`sampleID` as prefix
`bam` as bam file path
  * output:  
`SampleID.my.count.rds`
2. create main matrix of read count data of all sample
3. run CallCNVs for each sample
