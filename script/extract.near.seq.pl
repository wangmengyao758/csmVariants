#!/usr/bin/perl -w

my $hg38="/home/wangmengyao/GATK/database/Homo_sapiens_assembly38.fasta";
my $samtools="/apps/software/SAMtools/1.15.1-GCC-11.2.0/bin/samtools";
open OUT, ">extract.near.dbsnp.seq.fa";
open IN, "summary.frequency.txt" or die "$!\n";
while(<IN>){
	chomp;
	my ($gene,$chr,$pos) = (split)[0,1,2];
	my $s=$pos - 10;
	my $e=$pos + 10;
	my $re = "chr"."$chr".":"."$s"."-"."$e";
	my $id = "$gene"."_"."$chr"."_"."$pos";
	my $seq = `$samtools faidx $hg38 $re |grep  -v ">"`;
	chomp $seq;
	print OUT ">$id\n$seq\n";
}
close IN;
close OUT;

