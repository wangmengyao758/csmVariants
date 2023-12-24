#!/usr/bin/perl -w

my %hash;
open IN, "TCGA.merge.somatic.snv.dbsnp.regulatory_features.bed" or die "$!\n";
while(<IN>){
	chomp;
	my ($chr,$pos,$type) = (split)[0,1,-2];
	my $key = "$chr\t$pos";
	next if($type eq ".");
	$hash{$key} = $type;
}
close IN;

open CI, "summary.tcga.cancertype.dbsnp.table.info.txt" or die "$!\n";
open OUT, ">summary.tcga.cancertype.dbsnp.table.regulatory.info.txt";
my $head = <CI>;
chomp $head;
print OUT "$head\tRegulatory\n";
while(<CI>){
	chomp;
	my ($chr,$pos) = (split)[2,3];
	my $key = "chr$chr\t$pos";
	my $reg = "none";
	if(exists $hash{$key}){$reg = $hash{$key}}
	print OUT "$_\t$reg\n";
}
close CI;
close OUT;
