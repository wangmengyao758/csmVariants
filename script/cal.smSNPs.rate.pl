#!/usr/bin/perl -w

my (%hashs,%hash1,%hash2,%hash3,%hash);

open IN, "/home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.maf" or die "$!\n";
while(<IN>){
        chomp;
	next if(/^#/);
	next if(/^Hugo_Symbol/);
        my ($chr,$pos,$type,$ref,$alt,$sample) = (split /\t/,$_)[4,5,9,10,12,15];
        next if($type ne "SNP");
	my $sa = substr($sample,0,12);
        $hash1{$sa} += 1;
}
close IN;

open IN1, "smSNPs.context.txt" or die "$!\n";
while(<IN1>){
	chomp;
	my ($id,$context) = (split)[0,-1];
	my ($i1,$i2) = (split /_/,$id)[0,1];
	my ($chr,$pos) = (split /:/,$i1)[0,1];
	my ($ref,$alt) = (split /\//,$i2)[0,1];
	my $key = "$chr\t$pos";
	$hash{$key} = "$context\t$ref\t$alt";
}
close IN1;
my $s=0;
open IN3, "/home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.dbsnp.maf" or die "$!\n";
<IN3>;
while(<IN3>){
        chomp;
        my ($chr,$pos,$type,$ref,$alt,$sample) = (split /\t/,$_)[4,5,9,10,12,15];
        next if($type ne "SNP");
	my $sa = substr($sample,0,12);
	my $key = "$chr\t$pos";
	#my $context = (split /\t/,$hash{$key})[0];
	$hash2{$sa} += 1;
	my $tt = 0;
	if(exists $hash{$key}){
		my $context = (split /\t/,$hash{$key})[0];
		if($context =~ /CG/ || $context =~ /GC/){$tt = 1}
		if($ref eq "C" && $alt eq "T" && $tt ==1){$hash3{$sa} += 1}
	}else{$s+=1}
}
close IN3;
print "$s\n";
open CI, "../tcga.clin.info.txt" or die "$!\n";
<CI>;
while(<CI>){
        chomp;
	s/\#N\/A/NA/g;
	my ($id,$type,$age,$sex,$stage,$os,$ostime,$pfs,$pfstime) = (split /\t/,$_)[1,2,3,4,6,25,26,31,32];
	my $Stage = "NA";
	$stage =~ s/Stage\s//;
	$stage =~ s/[A|B|C]$//;
	if($stage =~ /\[/){$stage = "NA"}
	if($stage =~ /NOS/){$stage = "NA"}
	if($stage =~ /IS/){$stage = "NA"}
	if($stage =~ /0/){$stage = "NA"}
	if($stage eq "X"){$stage = "IV"}
	#if($stage eq "I" || $stage eq "IA" || $stage eq "IB"){$Stage = "I"}
	#if($stage eq "II" || $stage eq "IIA" || $stage eq "IIB" || $stage eq "IIC"){$stage = "II"}
	
	$Stage = $stage;
        $hashs{$id} = "$type\t$age\t$sex\t$Stage\t$os\t$ostime\t$pfs\t$pfstime";

	#my ($id,$type,$age,$sex) = (split)[1,2,3,4];
	#$hashs{$id} = "$type\t$age\t$sex";
}
close CI;

open OUT, ">tcga.smSNPs.rate.txt";
print OUT "Sample\tsmSNPs\ttotal\tsmSNPs_rate\tCpG.CtoT\tCpG.CtoT.rate\tcancer\tage\tgender\tStage\tOS\tOStime\tPFS\tPFStime\n";
foreach my $sa(sort keys %hash1){
     if(exists $hashs{$sa}){
	     my $a = $hash1{$sa};
	     my ($b,$c) = (0,0) ;
	     if(exists $hash2{$sa}){$b = $hash2{$sa};}
	     if(exists $hash3{$sa}){$c = $hash3{$sa};}
	     next if($a == 0);
             my $rate = $b / $a;
	     my $rate2 = $c / $a;
	     print OUT "$sa\t$b\t$a\t$rate\t$c\t$rate2\t$hashs{$sa}\n";
     }
}
close OUT;
