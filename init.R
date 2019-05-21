args=commandArgs(T)
#ref='hg19.fa'
ref=args[1]
library(ExomeDepth)
data(exons.hg19)
data(exons.hg19.X)
source("getBins.R")
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
