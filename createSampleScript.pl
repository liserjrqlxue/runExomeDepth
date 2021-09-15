#!/bin/env perl
#
use strict;
use warnings;
die "usage:perl $0 <input.list> <workdir> <outdir>\n",unless @ARGV==3;
my$list=$ARGV[0] or die"$0 input.list\n";

print "#!/bin/bash\nexport PATH=/share/backup/wangyaoshen/local/bin:\$PATH\n";
my@samples;
open IN,"< $list" or die$!;
#open OUT,"> $ARGV[2]" or die$!;
open OUT,"> $ARGV[2]/sample.list.checked" or die$!;
my %gender;
while(<IN>){
  chomp;
  /^main_sample_num/ and next;
#  my($sampleID,$bam)=split /\t/,$_;
	my ($sampleID)=split /\t/,$_;
	my $bam="$ARGV[1]/$sampleID/bwa/$sampleID.bqsr.bam";
	open GD,"$ARGV[1]/$sampleID/coverage/gender.txt",or die $!;
	my $gender;
	while (<GD>){
		my @tmp=split/\s+/,$_;
		if ($tmp[7] eq "Male"){
			$gender="M";
		}elsif ($tmp[7] eq "Female"){
			$gender="F";
		}
	}
	unless (exists $gender{$sampleID}){
  if(-e $bam){
#    print OUT join("\t",$sampleID,$bam),"\n";
		print OUT join("\t",$sampleID,$gender),"\n";
    push @samples,$sampleID;
    print "Rscript run.getBamCount.R $sampleID $bam A $ARGV[1]/$sampleID/cnv/ &\n";
		print "Rscript run.getBamCount.R $sampleID $bam $gender $ARGV[1]/$sampleID/cnv/ &\n";
#    print "Rscript run.getBamCount.X.R $sampleID $bam &\n";
  }else{
    print STDERR "skip $sampleID : can not find $bam\n";
  }
	}
	$gender{$sampleID}=$gender;
}
close IN;
close OUT;
print "wait\n";
print "Rscript run.getAllCounts.R A $ARGV[2]\n";
for(@samples){
	print "Rscript run.getCNVs.R $_ A $ARGV[2] &\n";
}
print "wait\n";
print "Rscript run.getAllCounts.R F $ARGV[2]\n";
print "Rscript run.getAllCounts.R M $ARGV[2]\n";
for(@samples){
  print "Rscript run.getCNVs.R $_ $gender{$_} &\n";
}
print "wait\n";
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
