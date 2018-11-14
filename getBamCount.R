getBamCount <- function(
    target,
    bam, 
    index = bam,
    min.mapq = 20, 
    read.width = 300, 
    prefix = 'exons.hg19'
    ) {
    
    my.param <- Rsamtools::ScanBamParam(
        flag = Rsamtools::scanBamFlag(
            isDuplicate = FALSE, 
            isPaired = TRUE, 
            isProperPair = TRUE, 
            isSecondaryAlignment = FALSE
            ),
        what = c("mapq", "pos", "isize"),
        )
    target =readRDS(paste0(prefix,".target.rds"))
    
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
