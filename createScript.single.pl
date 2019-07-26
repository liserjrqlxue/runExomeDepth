#!/bin/env perl
#
use strict;
use warnings;

use FindBin qw($Bin);

my$project="B2C_SGD";


my$sampleID=shift or
die"$0 sampleID gender bam outdir control\n";
my$gender=shift or
die"$0 sampleID gender bam outdir control\n";
my$bam=shift or
die"$0 sampleID gender bam outdir control\n";
my$outdir=shift or die"$0 sampleID gender bam outdir control\n";
my$control=shift or die"$0 sampleID gender bam outdir control\n";
system("mkdir -p $outdir/$sampleID")==0 or die$!;

open SL,"> $outdir/$sampleID/sample.list.checked" or die$!;
print SL join("\t",$sampleID,$gender),"\n";
close SL;

open SH,"> $outdir/$sampleID/run.sh" or die$!;
open LST,"> $outdir/$sampleID/all.CNV.calls.list" or die$!;



print SH
  "#!/bin/bash\n",
  "export LOCAL=/home/wangyaoshen/local\n",
  "export GCC=\$LOCAL/gcc-8.2.0\n",
  "export PATH=\$GCC/bin:\$LOCAL/bin:\$PATH\n",
  "export CPATH=\$GCC/include:\$LOCAL/include\n",
  "export LIBRARY_PATH=\$GCC/lib64:\$GCC/lib:\$LOCAL/lib64:\$LOCAL/lib:\$LIBRARY_PATH\n",
  "export LD_LIBRARY_PATH=\$GCC/lib64:\$GCC/lib:\$LOCAL/lib64:\$LOCAL/lib:\$LD_LIBRARY_PATH\n",
  "Bin=$Bin\n",
  "outdir=$outdir/$sampleID\n";

print SH "\n# reads count for each sample\n";
my$thread=0;
print SH "Rscript \$Bin/run.getBamCount.R $sampleID $bam A \$outdir \$Bin\n";
print SH "# call CNVs for each sample\n";
print SH "Rscript \$Bin/run.getCNVsFromControl.R $sampleID A \$outdir $control.A.my.count.rds\n";
print SH "\n# reads count for each sample with gender $gender\n";
print SH " Rscript \$Bin/run.getBamCount.R $sampleID $bam $gender \$outdir \$Bin\n";
print SH "# call CNVs for each sample\n";
print SH "Rscript \$Bin/run.getCNVsFromControl.R $sampleID $gender \$outdir $control.$gender.my.count.rds\n";
my$CNV_anno="/share/backup/wangyaoshen/src/CNV_anno";

print LST "$sampleID.A.CNV.calls.tsv\n";
print LST "$sampleID.$gender.CNV.calls.tsv\n";
close LST;

print SH
"# anno cnv\n",
"CNV_anno=$CNV_anno\n",
"perl $CNV_anno/script/add_cn_split_gene.batch.pl ",
"\$outdir/all.CNV.calls.list ",
"\$outdir/sample.list.checked ",
"\$CNV_anno/database/database.gene.list.NM ",
"\$CNV_anno/database/gene_exon.bed ",
"\$CNV_anno/database/OMIM/OMIM.xls ",
"\$outdir/all.CNV.calls.anno.withoutHGMD\n";
print SH
"perl $CNV_anno/script/add_HGMD_gross.pl ",
"\$outdir/all.CNV.calls.anno.withoutHGMD ",
"$CNV_anno/database/hgmd-gross_all-ex1-20190426.tsv ",
"\$outdir/all.CNV.calls.anno\n";

close SH;
