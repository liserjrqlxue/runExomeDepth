#!/bin/env perl
#
use strict;
use warnings;

my$project="B2C_SGD";


my$indir=shift or die"$0 indir outdir\n";
my$outdir=shift or die"$0 indir outdir\n";
my$tag=shift;
system("mkdir -p $outdir")==0 or die$!;

open IN,"< $indir/name.hash" or die"no $indir/name.hash:$!\n";
open SH,"> $outdir/run.sh" or die$!;
open OUT,"> $outdir/sample.list.checked" or die$!;
open LST,"> $outdir/all.CNV.calls.list" or die$!;

my@samples;
my%gender;
my%bamPath;

print SH "#!/bin/bash\nexport PATH=/share/backup/wangyaoshen/local/bin:\$PATH\n";
print SH "\n# reads count for each sample\n";
my$sampleN=0;
while(<IN>){
  chomp;
  my($sampleID,$gender)=split /\t/,$_;
  my$bam="$indir/$sampleID/bwa/$sampleID.final.bam";
  $bamPath{$sampleID}=$bam;
  if(-e $bam){
    $sampleN++;
    print OUT join("\t",$sampleID,$bam,$gender),"\n";
    push @samples,$sampleID;
    push @{$gender{$gender}},$sampleID;
    print SH "  Rscript run.getBamCount.R $sampleID $bam $outdir &\n";
  }else{
    print STDERR "# skip $sampleID : can not find $bam\n";
  }
}
print STDERR "# load $sampleN samples\n";
close IN;
close OUT;
print SH "wait\n";

print SH "Rscript run.getAllCounts.R $outdir/sample.list.checked $outdir\n";
print SH "# call CNVs for each sample\n";
for(@samples){
  print SH "  Rscript run.getCNVs.R $_ $outdir &\n";
  print LST "$_.CNV.calls.tsv\n";
}
print SH "wait\n";


for my$gender(keys%gender){
  print SH "\n# reads count for each sample with gender $gender\n";
  my@samples=@{$gender{$gender}};

  for my$sampleID(@samples){
    my$bam=$bamPath{$sampleID};
    print SH "  Rscript run.getBamCount.X.R $sampleID $bam $gender $outdir &\n";
  }
  print SH "wait\n";

  print SH "Rscript run.getAllCounts.X.R $outdir/sample.list.checked $gender $outdir\n";
  print SH "# call CNVs for each sample\n";
  for(@samples){
    print SH "  Rscript run.getCNVs.X.R $_ $gender $outdir &\n";
    print LST "$_.$gender.CNV.calls.tsv\n";
  }
  print SH "wait\n";
}
close LST;

my$CNV_anno="/share/backup/wangyaoshen/src/CNV_anno";
print SH
"# anno cnv\n",
"perl $CNV_anno/script/add_cn_split_gene.batch.pl ",
"$outdir/all.CNV.calls.list ",
"$outdir/sample.list.checked ",
"$CNV_anno/database/database.gene.list.NM ",
"$CNV_anno/database/gene_exon.bed ",
"$CNV_anno/database/OMIM/OMIM.xls ",
"$outdir/all.CNV.calls.anno\n";

print STDERR "# submit cmd:\nqsub -cwd -l vf=".($sampleN*2)."G,p=$sampleN -P $project -N ExomeDepth.$tag $outdir/run.sh\n";
close SH;
__END__
#样品编号	样品比对结果bam文件路径
15D6652318	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/15D6652318/bwa/15D6652318.final.bam
16D0144787-A	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D0144787-A/bwa/16D0144787-A.final.bam
16D1725954	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1725954/bwa/16D1725954.final.bam
16D1725967	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1725967/bwa/16D1725967.final.bam
16D1730076	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/16D1730076/bwa/16D1730076.final.bam
17D0005933	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0005933/bwa/17D0005933.final.bam
17D0087302	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0087302/bwa/17D0087302.final.bam
17D0120297	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0120297/bwa/17D0120297.final.bam
17D0120304	/ifs7/B2C_SGD/PROJECT/PP12_Project/WES/20180713_all/17D0120304/bwa/17D0120304.final.bam
