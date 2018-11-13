# ExomeDepth for WES CNV calling

## simple run
```bash
Rscript run.R bam.list bai.list
```

## TO-DO:
### speed up
1. run getBamCounts parallel for each sample
2. create main matrix of read count data of all sample
3. run CallCNVs for each sample
