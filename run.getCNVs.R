args=commandArgs(T)
sampleID=args[1]
tag=args[2]
outdir=args[3]

suffix=".my.count.rds"
suffix=paste0(".",tag,suffix)

library(ExomeDepth)
if(tag=="A"){
	data(exons.hg19)
	exons=exons.hg19
}else{
	data(exons.hg19.X)
	exons=exons.hg19.X
}
data(Conrad.hg19)
exons.hg19.GRanges <- GenomicRanges::GRanges(
    seqnames = exons$chromosome,
    IRanges::IRanges(start=exons$start,end=exons$end),
    names = exons$name
    )
my.counts.matrix <- readRDS(file=paste0(outdir,"/",'all',suffix))
samples <- colnames(my.counts.matrix)
for(i in 1:length(samples)){
    if(sampleID!=samples[i]){
        next
    }
    my.choice <- select.reference.set(
        test.counts = my.counts.matrix[,i],
        reference.counts = my.counts.matrix[,-i],
        bin.length = (exons$end - exons$start)/1000,
        #n.bins.reduced = 10000
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
        chromosome = exons$chromosome,
        start = exons$start+1,
        end = exons$end,
        name = exons$name
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
    output.file <- paste0(outdir,"/",sampleID,".",tag,'.CNV.calls.tsv')
    write.table(file=output.file,x=all.exons@CNV.calls,row.names=F,sep="\t",quote=F)
    saveRDS(all.exons,file=paste0(outdir,"/",sampleID,".",tag,'.all.exons.rds'))
}
