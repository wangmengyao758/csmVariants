#!/usr/bin/perl -w
my %hash;
open CI, "../summary.frequency.txt" or die "$!\n";
while(<CI>){
	chomp;
	my ($gene,$chr,$pos,$count,$type) = (split)[0,1,2,3,6];
	my $key = "$gene"."_"."$chr"."_"."$pos";
	$hash{$key} = "$count\t$type";
}
close CI;

open CI2, "/home/wangmengyao/data3/dbsnp/legoplot/smSNPs.context.txt" or die "$!\n";
<CI2>;
while(<CI2>){
	chomp;
	my ($id,$context) = (split)[0,-1];
	my ($key,$mut) = (split /_/,$id)[0,1];
	$key =~ s/:/_/;
	my ($ref,$alt) = (split /\//,$mut)[0,1];
	$hashs{$key} = "$ref\t$alt\t$context";
}
close CI2;

open IN, "smSNPs.Palindromes.result.filter.txt" or die "$!\n";
my $head = <IN>; chomp $head;
open OUT, ">smSNPs.Palindromes.result.filter.anno.txt";
print OUT "$head\tcount\ttype\tref\talt\tcontext\n";
while(<IN>){
	chomp;
	my $id = (split)[0];
	my ($gene,$chr,$pos) = (split /_/,$id)[0,1,2];
	my $kk = "chr$chr"."_"."$pos";
	my ($count,$type,$ref,$alt,$context) = ("-","-","-","-","-");
	if(exists $hash{$id}){($count,$type) = (split /\t/,$hash{$id})[0,1]}
	if(exists $hashs{$kk}){($ref,$alt,$context) = (split /\t/,$hashs{$kk})[0,1,2];}
	if($ref eq "A"){$ref = "T"; $alt =~ tr/ATCG/TAGC/; $context =~ tr/ATCG/TAGC/; $context = reverse $context;}
	if($ref eq "G"){$ref = "C"; $alt =~ tr/ATCG/TAGC/; $context =~ tr/ATCG/TAGC/; $context = reverse $context;}

	print OUT "$_\t$count\t$type\t$ref\t$alt\t$context\n";
}
close IN;
close OUT;

