#!/usr/bin/perl -w
my %hash;


open IN1, "summary.tcga.cancertype.dbsnp.table.info.txt" or die "$!\n";
#open OUT, ">summary.tcga.cancertype.dbsnp.table.info.muttype.txt";
open OUT1, ">TCGA.smSNPs.for.CMplot.txt";
print OUT1 "CHR	SNP	BP	A1	TEST	NMISS	OR	STAT	Age	Gender\tGene\tType\n";
my $header = <IN1>;
chomp $header;
while(<IN1>){
	chomp;
	my @temps = split;
        my $chr = $temps[2];
	next if($chr eq "X");
	next if($chr eq "Y");
	my $pos = $temps[3];
	my $gene = $temps[1];
	my $count = $temps[4];
	my $age = $temps[69];
	my $sex = $temps[68];
	my $type = $temps[11];
	my $id = "chr"."$chr".":".$pos;
	print OUT1 "$chr\t$id\t$pos\tNA\tADD\t$count\tNA\tNA\t$age\t$sex\t$gene\t$type\n";
}
close IN1;
close OUT1;

