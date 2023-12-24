#!/usr/bin/perl -w

my (%hash,%hashs,%hashp);
open IN, "summary.tcga.cancertype.dbsnp.table.fisher.test.txt" or die "$!\n";
<IN>;
while(<IN>){
	chomp;
	my ($gene,$chr,$pos,$p_sex,$p_age) = (split)[0,1,2,14,15];
	my $key = "$chr\t$pos";
	$hash{$key} = "$p_sex\t$p_age";
}
close IN;

open (IN1, "gunzip -c /home/wangmengyao/data3/dbsnp/hg38/common_all_20180418.vcf.gz|") or die "$!\n";
while(<IN1>){
	chomp;
	next if(/^#/);
	my ($chr,$pos,$snp,$ref,$alt) = (split)[0,1,2,3,4];
	$chr =~ s/chr//;
	my $key = "$chr\t$pos";
	$hashs{$key} = $snp;
}
close IN1;

open IN2, "TCGA.merge.somatic.snv.dbsnp.phastCons100way.UCSC.annot.txt" or die "$!\n";
<IN2>;
while(<IN2>){
	chomp;
	my ($chr,$pos,$score) = (split)[0,1,-1];
        $chr =~ s/chr//;
	my $key = "$chr\t$pos";
	$hashp{$key} = $score;
}
close IN2;

open FI, "summary.tcga.cancertype.dbsnp.table.txt" or die "$!\n";
open OUT, ">summary.tcga.cancertype.dbsnp.table.info.txt";
my $head = <FI>;
chomp $head;
print OUT "SNPid\t$head\tSex.p\tAge.p\tphastCons.score\n";
while(<FI>){
	chomp;
        my ($chr,$pos) = (split)[1,2];
	my $key = "$chr\t$pos";
	my ($snp,$p1,$p2,$score) = ("NA","NA","NA","NA");
	if(exists $hashs{$key}){$snp = $hashs{$key}}
	if(exists $hash{$key}){($p1,$p2) = (split /\t/,$hash{$key})[0,1]}
	if(exists $hashp{$key}){$score = $hashp{$key}}
	print OUT "$snp\t$_\t$p1\t$p2\t$score\n";
}
close FI;
close OUT;
