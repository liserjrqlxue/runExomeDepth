args=commandArgs(T)
sampleID=args[1]
outdir=args[2]
library(ExomeDepth)
data(exons.hg19)
data(Conrad.hg19)
exons.hg19.GRanges <- GenomicRanges::GRanges(
    seqnames = exons.hg19$chromosome,
    IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
    names = exons.hg19$name
    )
my.counts.matrix <- readRDS(file=paste0(outdir,"/",'all',".my.counts.rds"))
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
        column.name = 'exons.hg19'
        )
    output.file <- paste0(outdir,"/",sampleID,'.CNV.calls.tsv')
    write.table(file=output.file,x=all.exons@CNV.calls,row.names=F,sep="\t",quote=F)
    saveRDS(all.exons,file=paste0(outdir,"/",sampleID,'.all.exons.rds'))
}
