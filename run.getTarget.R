library(ExomeDepth)
data(exons.hg19)
source("getBins.R")
getBins(
    bed.frame = exons.hg19,
    include.chr = T
)
