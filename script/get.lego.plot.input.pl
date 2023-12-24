#!/usr/bin/perl -w
my $head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT";
my (%hash1,%hash2,%hash3);
open (IN, "gunzip -c /home/wangmengyao/data3/dbsnp/hg38/common_all_20180418.vcf.gz|" ) or die "$!\n";
open OUT1, ">dbsnps.context.hg38.vcf";
print OUT1 "$head\n";
while(<IN>){
        chomp;
        next if(/^#/);
        my ($chr,$pos,$snp,$ref,$alt) = (split /\t/,$_)[0,1,2,3,4];
        my $a = length($ref);
        my $b = length($alt);
        next if($a!=1);
        next if($b!=1);
        $hash1{$snp} = "$ref\t$alt";
	print OUT1 "chr$chr\t$pos\t$snp\t$ref\t$alt\n";
}
close IN;
close OUT1;

open OUT2, ">tcga.context.hg38.vcf";
print OUT2 "$head\n";
open (IN2, "gunzip -c /home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.maf.gz|") or die "$!\n";
<IN2>;
while(<IN2>){
	chomp;
	my ($chr,$pos,$type,$ref,$alt) = (split /\t/,$_)[4,5,9,10,12];
	next if($type ne "SNP");
	my $key = "$chr\t$pos";
	next if(exists $hash2{$key});
	$hash2{$key} = 1;
	print OUT2 "$chr\t$pos\t.\t$ref\t$alt\n";
}
close IN2;
close OUT2;

open OUT3, ">smsnps.context.hg38.vcf";
print OUT3 "$head\n";
open IN3, "/home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.dbsnp.maf" or die "$!\n";
<IN3>;
while(<IN3>){
	chomp;
	my ($chr,$pos,$type,$ref,$alt) = (split /\t/,$_)[4,5,9,10,12];
	next if($type ne "SNP");
	my $key = "$chr\t$pos";
	next if(exists $hash3{$key});
	$hash3{$key} = 1;
	print OUT3 "$chr\t$pos\t.\t$ref\t$alt\n";
}
close IN3;
close OUT3;
