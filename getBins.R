getBins <- function(
    bed.frame = NULL,
    bed.file = NULL,
    read.width = 300, 
    include.chr = FALSE, 
    referenceFasta = NULL,
    prefix = 'exons.hg19',
    ) {

    if (is.null(bed.frame)) {
        if (is.null(bed.file)) {
            stop("If no bed data frame is provided there must be a link to a bed file")
        }
        bed.frame <- read.delim(file = bed.file, header =  FALSE, stringsAsFactors = FALSE)
    }

    names(bed.frame)[1] <- 'seqnames'
    names(bed.frame)[2] <- 'start'
    names(bed.frame)[3] <- 'end'

    if (include.chr) {
        if (sum(grepl(pattern = '^chr', bed.frame$seqnames) > 0)) {
            warning('The option include.chr == TRUE adds the chr prefix to the chromosome name but it looks like the chromosome names already have a chr prefix. The argument to getBamCounts is probably an error.')
        }
        bed.frame$seqnames <- paste('chr', bed.frame$seqnames, sep = '')
    }

    chr.names.used <- unique(as.character(bed.frame$seqnames))
    chr.levels <- c(as.character(seq(1, 22)), subset( chr.names.used, ! chr.names.used %in% as.character(seq(1, 22))))

    ####specifying the levels is important here to not mess up the order
    bed.frame$seqnames <- factor(bed.frame$seqnames, levels = chr.levels)  

    ##order the data frame by position
    bed.frame <- bed.frame[ order(bed.frame$seqnames, bed.frame$start + bed.frame$end), ]  

    target <- GenomicRanges::GRanges(
        seqnames = bed.frame$seqnames,
        IRanges::IRanges(start=bed.frame$start+1,end=bed.frame$end)
        )
    saveRDS(target,file=paste0(prefix,".target.rds"))

    rdata <- IRanges::RangedData(
        space = GenomicRanges::seqnames(target),
        ranges= GenomicRanges::ranges(target)
        )

    if ((ncol(bed.frame) >= 4) && (class(bed.frame[,4]) %in% c('character', 'factor'))) {
        ##add exon names if available
        row.names(rdata) <- make.unique(as.character(bed.frame[,4]))  
    }

    ############################################################################# add GC content
    if (!is.null(referenceFasta)) {
        message('Reference fasta file provided so ExomeDepth will compute the GC content in each window')
        target.dnastringset <- Rsamtools::scanFa(referenceFasta, target)

        getGCcontent <- function(x) {
            GC.count  <- Biostrings::letterFrequency(x,"GC")
            all.count <- Biostrings::letterFrequency(x,"ATGC")
            as.vector(ifelse(all.count==0,NA,GC.count/all.count))
        }
        rdata[["GC"]] <- getGCcontent(target.dnastringset)
    }
    saveRDS(rdata,file=paste0(prefix,".rdata.rds"))
    return(rdata)
}
