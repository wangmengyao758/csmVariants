#!/usr/bin/perl -w

my %hash;

open (IN, "gunzip -c hg38/common_all_20180418.vcf.gz|" ) or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/);
	my ($chr,$pos,$ref,$alt) = (split)[0,1,3,4];
	$chr =~ s/chr//;
	my $key = "$chr\t$pos";
	$hash{$key} = "$ref\t$alt";
}
close IN;
#open (IN2, "gunzip -c hg38/mbiobank_ChinaMAP.phase1.vcf.gz |") or die "$!\n";
#while(<IN2>){
	#chomp;
	#next if(/^#/);
	#my ($chr,$pos,$ref,$alt) = (split)[0,1,3,4];
	#$chr =~ s/chr//;
	#my $key = "$chr\t$pos";
	#$hash{$key} = "$ref\t$alt";
	#}
#close IN2;
my $head;
open OUT, ">TCGA.merge.somatic.snv.dbsnp.maf";
open CI, "TCGA.merge.somatic.snv.maf" or die "$!\n";
while(<CI>){
	chomp;
	next if(/^#/);
	if(/^Hugo_Symbol/){
		if(!$head){$head=$_;print OUT "$head\n"}
	}else{
	  my ($chr,$pos,$ref,$alt,$tref,$talt,$nref,$nalt) = (split /\t/,$_)[4,5,10,12,40,41,43,44];
	  $chr =~ s/chr//;
	  next if($nalt >=3);
	  my $tsum = $talt + $tref;
	  next if($tsum < 5);
	  my $tfre = $talt / $tsum;
	  next if($tfre < 0.05);
	  my $key = "$chr\t$pos";
	  my $aa = "$ref"."/"."$alt";
	  if(exists $hash{$key}){print OUT "$_\n"}
        }
}
close CI;
close OUT;

