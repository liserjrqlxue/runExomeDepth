getBamCount <- function(
    target,
    bam, 
    index = bam,
    min.mapq = 20, 
    read.width = 300
    ) {
    
    t <- countBamInGRanges.exomeDepth(
        bam.file = bam, 
        index = index, 
        granges = target, 
        min.mapq = min.mapq, 
        read.width = read.width
        )
    message("Number of counted fragments : ", sum(t))
    return(t)
}
