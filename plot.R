args=commandArgs(T)
library(ExomeDepth)
sampleID=args[1]
rds=args[2]
gene=args[3]
chr=args[4]
start=as.numeric(args[5])
stop=as.numeric(args[6])

regLen=stop-start
regLen
d=2
xlim=c(start-regLen/d,stop+regLen/d)
xlim
pdf(paste(sampleID,gene,"pdf",sep="."),width=16,height=9)
all.exons=readRDS(rds)
plot(
  all.exons,
  sequence=chr,
  xlim=xlim,
  count.threshold=10,
  main=paste(sampleID,gene,sep=" "),
  cex.lab=0.8,
  with.gene=TRUE
)
dev.off()
message(paste0("http://192.168.3.4:9091/ifs7/B2C_SGD/PROJECT/PP12_Project/wangyaoshen/ExomeDepth/",paste(sampleID,gene,"pdf",sep=".")))
