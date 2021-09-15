#!/bin/env perl
#
use strict;
use warnings;

use FindBin qw($Bin);

#my$project="B2C_SGD";

#my$MaxThread=shift or die "$0 maxThread indir input outdir\n";;
my $indir  = shift or die "$0 indir input outdir\n";
my $input  = shift or die "$0 indir input outdir\n";
my $outdir = shift or die "$0 indir input outdir\n";
my $tag    = shift;
system("mkdir -p $outdir") == 0 or die $!;

open IN,  "< $input"                      or die "no $input:$!\n";
open SH,  "> $outdir/run.sh"              or die $!;
open OUT, "> $outdir/sample.list.checked" or die $!;
open ML,  "> $outdir/M.sample.list"       or die $!;
open FL,  "> $outdir/F.sample.list"       or die $!;
open LST, "> $outdir/all.CNV.calls.list"  or die $!;

my @samples;
my %gender;
my %bamPath;
my %hash;

print SH
  "#!/bin/bash\n",
  "set -exuo pipefail\n",

#"export LOCAL=/home/wangyaoshen/local\n",
#"export GCC=\$LOCAL/gcc-8.2.0\n",
#"export PATH=\$GCC/bin:\$LOCAL/bin:\$PATH\n",
#"export CPATH=\$GCC/include:\$LOCAL/include\n",
#"export LIBRARY_PATH=\$GCC/lib64:\$GCC/lib:\$LOCAL/lib64:\$LOCAL/lib:\$LIBRARY_PATH\n",
#"export LD_LIBRARY_PATH=\$GCC/lib64:\$GCC/lib:\$LOCAL/lib64:\$LOCAL/lib:\$LD_LIBRARY_PATH\n",
  "Bin=$Bin\n", "outdir=$outdir\n";

print SH "\n# reads count for each sample\n";
my $sampleN = 0;
my $thread  = 0;
while (<IN>) {
    chomp;
    next if ( $_ =~ /main_sample_num/ );

    #  my($sampleID,$gender)=split /\t/,$_;
    my ($sampleID) = split /\t/, $_;

    #  my$bam="$indir/$sampleID/bwa/$sampleID.final.bam";
    my $bam = "$indir/$sampleID/bwa/$sampleID.bqsr.bam";
    $bamPath{$sampleID} = $bam;
    unless ( exists $hash{$sampleID} ) {

        #if(-e $bam){
        $sampleN++;
        open F, "$indir/$sampleID/coverage/gender.txt", or die $!;
        my $line = <F>;
        chomp $line;
        my $gender;
        my @gender = split /\s+/, $line;
        if ( $gender[7] eq "Male" ) {
            $gender = "M";
            print ML join( "\t", $sampleID, $gender, $bam ), "\n";
        }
        if ( $gender[7] eq "Female" ) {
            $gender = "F";
            print FL join( "\t", $sampleID, $gender, $bam ), "\n";
        }
        $hash{$sampleID} = $gender;
        print OUT join( "\t", $sampleID, $gender, $bam ), "\n";
        push @samples, $sampleID;
        push @{ $gender{$gender} }, $sampleID;
        #print SH "  Rscript \$Bin/run.getBamCount.R $sampleID $bam A \$outdir \$Bin &\n";
    }
}
print STDERR "# load $sampleN samples\n";
close IN;
close OUT;
close ML;
close FL;
#print SH "wait\n";
print SH "cat $outdir/sample.list.checked|gargs -l $outdir/getBamCount.A.log -p 10 'Rscript $Bin/run.getBamCount.R {0} {2} A $outdir $Bin'\n";

print SH
  "Rscript \$Bin/run.getAllCounts.R \$outdir/sample.list.checked A \$outdir\n";
print SH "# call CNVs for each sample\n";
for (@samples) {
    #print SH "  Rscript \$Bin/run.getCNVs.R $_ A \$outdir&\n";
    print LST "$_.A.CNV.calls.tsv\n";
}
print SH "cat $outdir/sample.list.checked|gargs -l $outdir/getCNVs.A.log -p 10 'Rscript $Bin/run.getCNVs.R {0} A $outdir'\n";
#print SH "wait\n";

for my $gender (sort keys %gender ) {
    print SH "\n# reads count for each sample with gender $gender\n";
    my @samples = @{ $gender{$gender} };

    for my $sampleID (@samples) {
        my $bam = $bamPath{$sampleID};
        #print SH "  Rscript \$Bin/run.getBamCount.R $sampleID $bam $gender \$outdir \$Bin &\n";
    }
    #print SH "wait\n";
    print SH "cat $outdir/$gender.sample.list|gargs -l $outdir/getBamCount.$gender.log -p 10 'Rscript $Bin/run.getBamCount.R {0} {2} $gender $outdir $Bin'\n";

    print SH
"Rscript \$Bin/run.getAllCounts.R \$outdir/sample.list.checked $gender \$outdir\n";
    print SH "# call CNVs for each sample\n";
    for (@samples) {
        #print SH "  Rscript \$Bin/run.getCNVs.R $_ $gender \$outdir &\n";
        print LST "$_.$gender.CNV.calls.tsv\n";
    }
    #print SH "wait\n";
    print SH "cat $outdir/$gender.sample.list|gargs -l $outdir/getCNVs.$gender.log -p 10 'Rscript $Bin/run.getCNVs.R {0} $gender $outdir'\n";
}
close LST;

#my$CNV_anno="/share/backup/wangyaoshen/src/CNV_anno";
print "$Bin/../CNV_anno\n";
my $CNV_anno = "$Bin/../CNV_anno";
print SH
  "# anno cnv\n",
  "CNV_anno=$CNV_anno\n";
for (@samples) {
    #print SH "perl $CNV_anno/script/add_cn_split_gene.pl $_ $outdir/$_.A.CNV.calls.tsv,$outdir/$_.$hash{$_}.CNV.calls.tsv $hash{$_} $CNV_anno/database/database.gene.list.NM $CNV_anno/database/gene_exon.bed $CNV_anno/database/OMIM/OMIM.xls $outdir/$_.CNV.calls.anno.withoutHGMD\n";
    #print SH "perl $CNV_anno/script/add_HGMD_gross.pl $outdir/$_.CNV.calls.anno.withoutHGMD $CNV_anno/database/hgmd-gross_all-ex1-20190426.tsv $outdir/$_.CNV.calls.anno\n";
}
print SH "cat $outdir/sample.list.checked|gargs -l $outdir/add_cn_split_gene.log -p 10 'perl $CNV_anno/script/add_cn_split_gene.pl {0} $outdir/{0}.A.CNV.calls.tsv,$outdir/{0}.{1}.CNV.calls.tsv {1} $CNV_anno/database/database.gene.list.NM $CNV_anno/database/gene_exon.bed $CNV_anno/database/OMIM/OMIM.xls $outdir/{0}.CNV.calls.anno.withoutHGMD'\n";
print SH "cat $outdir/sample.list.checked|gargs -l $outdir/add_HGMD_gross.log -p 10 'perl $CNV_anno/script/add_HGMD_gross.pl $outdir/{0}.CNV.calls.anno.withoutHGMD $CNV_anno/database/hgmd-gross_all-ex1-20190426.tsv $outdir/{0}.CNV.calls.anno'\n";

=head
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
=cut

#print STDERR "# submit cmd:\nqsub -cwd -l vf=".($sampleN*2.2)."G,p=$sampleN -P $project -N ExomeDepth.$tag $outdir/run.sh\n";
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
