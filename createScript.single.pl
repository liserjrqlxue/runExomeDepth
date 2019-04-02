#!/bin/env perl
#
use strict;
use warnings;

use FindBin qw($Bin);

my$project="B2C_SGD";


my$sampleID=shift or die$!; 
my$gender=shift or die$!;
my$indir=shift or die"$0 indir outdir\n";
my$outdir=shift or die"$0 indir outdir\n";
my$tag=shift;
system("mkdir -p $outdir")==0 or die$!;

open IN,"< $indir/name.hash" or die"no $indir/name.hash:$!\n";
open SH,"> $outdir/$sampleID/run.sh" or die$!;
open LST,"> $outdir/$sampleID/all.CNV.calls.list" or die$!;

my%gender;
my%bamPath;

print SH 
  "#!/bin/bash\nexport PATH=/share/backup/wangyaoshen/local/bin:\$PATH\n",
  "Bin=$Bin\n",
  "outdir=$outdir/$sampleID\n";

print SH "\n# reads count for each sample\n";
my$thread=0;
my$bam="$indir/$sampleID/bwa/$sampleID.final.bam";
print SH "Rscript \$Bin/run.getBamCount.R $sampleID $bam A \$outdir \$Bin\n";
print SH "# call CNVs for each sample\n";
print SH "Rscript \$Bin/run.getCNVsFromControl.R $sampleID A \$outdir\n";
print SH "\n# reads count for each sample with gender $gender\n";
print SH " Rscript \$Bin/run.getBamCount.R $sampleID $bam $gender \$outdir \$Bin\n";
print SH "# call CNVs for each sample\n";
print SH "Rscript \$Bin/run.getCNVsFromControl.R $sampleID $gender \$outdir \n";
my$CNV_anno="/share/backup/wangyaoshen/src/CNV_anno";

print LST "$outdir/$sampleID/$sampleID.A.CNV.calls.tsv\n";
print LST "$outdir/$sampleID/$sampleID.$gender.CNV.calls.tsv\n";
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
"\$outdir/all.CNV.calls.anno\n";

close SH;
__END__
#样品编号  样品比对结果bam文件路径
15D6652318  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/15D6652318/bwa/15D6652318.final.bam
16D0144787-A  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D0144787-A/bwa/16D0144787-A.final.bam
16D1725954  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1725954/bwa/16D1725954.final.bam
16D1725967  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1725967/bwa/16D1725967.final.bam
16D1730076  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1730076/bwa/16D1730076.final.bam
17D0005933  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0005933/bwa/17D0005933.final.bam
17D0087302  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0087302/bwa/17D0087302.final.bam
17D0120297  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0120297/bwa/17D0120297.final.bam
17D0120304  /ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0120304/bwa/17D0120304.final.bam
