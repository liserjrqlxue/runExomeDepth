#!/share/backup/wangyaoshen/local/bin/R
args=commandArgs(T)
bamList=args[1]
baiList=args[2]

library(ExomeDepth)

bam.list<-read.table(bamList,stringsAsFactors=F)[,1]
bai.list<-read.table(baiList,stringsAsFactors=F)[,1]

fasta="hg19_chM_male_mask.fa"

# 3 Create count data from BAM files
## 3.1 Count for autosomal chromosomes

# target region
data(exons.hg19)

my.counts <- getBamCounts(
	bed.frame = exons.hg19,
	bam.files = bam.list,
	index.files = bai.list,
	include.chr = T,
	referenceFasta = fasta
)
saveRDS(my.counts,file=paste0(bamList,"my.counts.rds"))

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
my.counts.dafr$chromosome <- 
	gsub(	##remove the annoying chr letters
		as.character(my.counts.dafr$space),
		pattern = 'chr',
		replacement = ''
	)
#print(head(my.counts.dafr))


## 3.2 Counts for chromosome X
#data(exons.hg19.X)
#head(exons.hg19.X)

# 5 Build the most appropriate reference set
#my.test 	<- my.counts$"16D0144787-A.final.bam"
#my.ref.samples 	<- c('X15D6652318.final.bam','X16D1725954.final.bam','X16D1725967.final.bam')
#my.reference.set <- as.matrix(my.counts.dafr[, my.ref.samples])
#my.choice 	<- select.reference.set (
#	test.counts = my.test,
#	reference.counts = my.reference.set,
#	bin.length = (my.counts.dafr$end - my.counts.dafr$start)/1000,
#	n.bins.reduced = 10000 # speed up
#)
## Optimization of the choice of aggregate reference set, this process can take some time
## Number of selected bins: 10000
## Warning in aod::betabin(data = data.for.fit, formula = as.formula(formula), : The data set contains at least one line with weight = 0.
#print(my.choice[[1]])
## [1] "Exome2" "Exome1" "Exome3"

#my.matrix <- as.matrix( my.counts.dafr[, my.choice$reference.choice, drop = FALSE])
#my.reference.selected <- apply(
#	X = my.matrix,
#	MAR = 1,
#	FUN = sum
#)


# 6 CNV calling
#all.exons <- new(
#	'ExomeDepth',
#	test = my.test,
#	reference = my.reference.selected,
#	formula = 'cbind(test, reference) ~ 1'
#)
## Now fitting the beta-binomial model on a data frame with 26547 rows : this step can take a few minutes.
## Warning in aod::betabin(data = data.for.fit, formula = as.formula(formula), : 
## The data set contains at least one line with weight = 0.
## Now computing the likelihood for the different copy number states

#all.exons <- CallCNVs(
#	x = all.exons,
#	transition.probability = 10^-4,
#	#chromosome = my.counts.dafr$space,
#	chromosome = my.counts.dafr$chromosome,
#	start = my.counts.dafr$start,
#	end = my.counts.dafr$end,
#	name = my.counts.dafr$names
#)
## Correlation between reference and tests count is 0.9898
## To get meaningful result, this correlation should really be above 0.97. If this is not the case,
## consider the output of ExomeDepth as less reliable (i.e. most likely a high false positive rate)
## Number of calls for chromosome 1 : 25
#head(all.exons@CNV.calls)
## start.p end.p type nexons start end chromosome
## 1 25 27 deletion 3 89553 91106 1
## 2 62 66 deletion 5 461753 523834 1
## 3 100 103 duplication 4 743956 745551 1
## 4 575 576 deletion 2 1569583 1570002 1
## 5 587 591 deletion 5 1592941 1603069 1
## 6 2324 2327 deletion 4 12976452 12980570 1
## id BF reads.expected reads.observed reads.ratio
## 1 chr1:89553-91106 12.40 224 68 0.304
## 2 chr1:461753-523834 9.82 363 190 0.523
## 3 chr1:743956-745551 7.67 201 336 1.670
## 4 chr1:1569583-1570002 5.53 68 24 0.353
## 5 chr1:1592941-1603069 13.90 1136 434 0.382
## 6 chr1:12976452-12980570 12.10 780 342 0.438

#output.file <- 'exome_calls.csv'
#write.csv(
#	file = output.file,
#	x = all.exons@CNV.calls,
#	row.names = FALSE
#)

# 7 Ranking the CNV calls by confidence level
#head(all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),])

# 8 Better annotation of CNV calls
#data(Conrad.hg19)
#head(Conrad.hg19.common.CNVs)
## GRanges object with 6 ranges and 1 metadata column:
## seqnames ranges strand | names
## <Rle> <IRanges> <Rle> | <factor>
## [1] 1 [10499, 91591] * | CNVR1.1
## [2] 1 [10499, 177368] * | CNVR1.2
## [3] 1 [82705, 92162] * | CNVR1.5
## [4] 1 [85841, 91967] * | CNVR1.4
## [5] 1 [87433, 89163] * | CNVR1.6
## [6] 1 [87446, 109121] * | CNVR1.7
## -------
## seqinfo: 23 sequences from an unspecified genome; no seqlengths

#levels(GenomicRanges::seqnames(Conrad.hg19.common.CNVs))
## [1] "1" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "2" "20" "21"
## [15] "22" "3" "4" "5" "6" "7" "8" "9" "X"

#all.exons <- AnnotateExtra(
#	x = all.exons,
#	reference.annotation = Conrad.hg19.common.CNVs,
#	min.overlap = 0.5,
#	column.name = 'Conrad.hg19'
#)
#print(head(all.exons@CNV.calls))

#exons.hg19.GRanges <- GenomicRanges::GRanges(
#	seqnames = exons.hg19$chromosome,
#	IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
#	names = exons.hg19$name
#)
#all.exons <- AnnotateExtra(
#	x = all.exons,
#	reference.annotation = exons.hg19.GRanges,
#	min.overlap = 0.0001,
#	column.name = 'exons.hg19'
#)

# 9 Visual display
#plot (
#	all.exons,
#	sequence = '19',
#	xlim = c(36326608 - 1e4, 36326663 + 1e4),
#	count.threshold = 20,
#	main = 'RHD gene',
#	cex.lab = 0.8,
#	with.gene = TRUE
#)
## Plotting the gene data


# 10 How to loop over the multiple samples

#### get the annotation datasets to be used later
data(Conrad.hg19)
exons.hg19.GRanges <- GenomicRanges::GRanges(
	seqnames = exons.hg19$chromosome,
	IRanges::IRanges(start=exons.hg19$start,end=exons.hg19$end),
	names = exons.hg19$name
)
### prepare the main matrix of read count data
my.counts.mat <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern = '.*.final.bam')])
nsamples <- ncol(my.counts.mat)
### start looping over each sample
for (i in 1:nsamples) {
	prefix=basename(bam.list)[i]
	message(paste("call cnvs for",prefix))
#### Create the aggregate reference set for this sample
	my.choice <- select.reference.set (
		test.counts = my.counts.mat[,i],
		reference.counts = my.counts.mat[,-i],
		bin.length = (my.counts.dafr$end - my.counts.dafr$start)/1000,
		n.bins.reduced = 10000
	)
	my.reference.selected <- apply(
		X = my.counts.mat[, my.choice$reference.choice, drop = FALSE],
		MAR = 1,
		FUN = sum
	)
	message('Now creating the ExomeDepth object')
	all.exons <- new(
		'ExomeDepth',
		test = my.counts.mat[,i],
		reference = my.reference.selected,
		formula = 'cbind(test, reference) ~ 1'
	)
################ Now call the CNVs
	all.exons <- CallCNVs(
		x = all.exons,
		transition.probability = 10^-4,
		#chromosome = my.counts.dafr$space,
		chromosome = my.counts.dafr$chromosome,
		start = my.counts.dafr$start,
		end = my.counts.dafr$end,
		name = my.counts.dafr$names
	)
########################### Now annotate the ExomeDepth object
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
	output.file <- paste(prefix, '.CNV.calls.tsv', sep = '')
	write.table(file = output.file, x = all.exons@CNV.calls, row.names = FALSE,sep="\t",quote=F)
}

# 11 Additional functions

## 11.2 Counting everted reads
#data(genes.hg19)
#everted <- count.everted.reads (bed.frame = genes.hg19,
#bam.files = bam.files.list,
#min.mapq = 20,
#include.chr = TRUE)
