#!/usr/bin/perl -w

my ($cancer,$snp) = @ARGV;

my $k = 10;
my $cfile = "$cancer.clin.info.txt";
my $mfile = "$cancer.merge.somatic.snv.maf";
my $ofile1 = "$cancer.merge.somatic.snv.$snp.mut.maf";
my $ofile2 = "$cancer.merge.somatic.snv.$snp.wild.maf";

my %hash;
open CI, "$cfile" or die "$!\n";
my $head = <CI>; chomp $head;
my @arrs = (split /\t/,$head);
for(my $i=10;$i<=$#arrs;$i+=1){if($arrs[$i] eq $snp){$k = $i}}
while(<CI>){
	chomp;
	my ($sample,$class) = (split /\t/,$_)[1,$k];
	$hash{$sample} = $class;
}
close CI;

open OUT1, ">$ofile1";
open OUT2, ">$ofile2";
open IN, "$mfile" or die "$!\n";
$head = <IN>;
print OUT1 "$head"; print OUT2 "$head";
while(<IN>){
	chomp;
	my $sa = (split /\t/,$_)[15];
	$sa = substr($sa,0,12);
	#	print "$sa\n";
	next if(!exists $hash{$sa});
	if($hash{$sa} == 0){print OUT2 "$_\n";}
	if($hash{$sa} == 1){print OUT1 "$_\n";}
	
}
close IN;
close OUT1;
close OUT2;


