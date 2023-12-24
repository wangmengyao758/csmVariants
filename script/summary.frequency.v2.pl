#!/usr/bin/perl -w

my (%hash, %hasht, %hashf,%hashi,%hasha);
my $total = 11252671;
my $dbsnp_sum = 37302986;
my $sum = 866176;
my $nsample = 10574;

open CI, "tcga.clin.info.txt" or die "$!\n";
<CI>;
while(<CI>){
	chomp;
        my ($id,$type,$age,$sex) = (split)[1,2,3,4];
	$hashi{$id} = "$type\t$age\t$sex";
}
close CI;

open (SN, "gunzip -c hg38/common_all_20180418.vcf.gz|" ) or die "$!\n";
while(<SN>){
        chomp;
        next if(/^#/);
        my ($chr,$pos,$id,$ref,$alt,$info) = (split)[0,1,2,3,4,7];
        $chr =~ s/chr//;
        my $key = "$chr\t$pos";
	$info =~ /CAF=([\d|\.|,]+);/;
	my $fres = $1;
	my @arrs = (split /,/,$fres);
	my @alts = (split /,/,$alt);
	$info =~ /VC=(\w+)/;
	my $type = $1;
	my $fre = 1 - $arrs[0];
	$fre = sprintf "%.5f",$fre;
	my $out = "$type\t$ref\t$fre\t$id";
	$hashf{$key} .= "$out\n";
	
}
close SN;

open IN, "TCGA.merge.somatic.snv.dbsnp.maf" or die "$!\n";
my $head = <IN>;
while(<IN>){
	chomp;
	my ($gene,$chr,$pos,$type,$snv,$ref,$alt,$sample) = (split /\t/,$_)[0,4,5,8,9,11,12,15];
	$chr =~ s/chr//;
	my $key = "$gene\t$chr\t$pos";
	my $value = "$type".";"."$sample";
	$hash{$key} .= "$sample;";
	$hasht{$key} .= "$type;";
	my $skey  = "$sample\t$key";
	$hasha{$skey} = "$type\t$snv\t$ref\t$alt";
}
close IN;

open COUT, ">summary.frequency.table.txt";
open OUT, ">summary.frequency.txt";
foreach my $key(sort keys %hash){
	my @arrs = split /;/,$hash{$key};
	my @types = split /;/,$hasht{$key};
	my ($gene,$chr,$pos) = (split /\t/,$key)[0,1,2];
	my $kk = "$chr\t$pos";
	#my $fre = (split /\t/, $hashf{$kk})[2];
        my %hashe; my $type;
	foreach my $tt(@types){$hashe{$tt} += 1;}
	foreach my $te(sort keys %hashe){my $aa = "$te".":"."$hashe{$te}"; $type .= "$aa;";}
	my $count = $#arrs + 1;
	my $fre1 = sprintf "%.5f", $count / $nsample;
        my $fre = 0; my $snpid="-";
	my @temps = (split /\n/,$hashf{$kk});

	my @samples = (split /;/,$hash{$key});
	foreach my $sample(@samples){
		my $skey = "$sample\t$key";
		my $sa = substr($sample,0,12);
		my ($snv,$r1,$a1) = (split /\t/,$hasha{$skey})[1,2,3];
		foreach my $temp(@temps){
			my ($tt2,$r2,$f2,$id2) = (split /\t/,$temp);
			
			if($tt2 eq "SNV" && $snv =~ /NP/){$fre = $f2; $snpid=$id2}
			elsif($tt2 eq "DIV" && ($snv =~ /INS/) || ($snv =~ /DEL/)){$fre = $f2; $snpid=$id2}
			
		}
		print COUT "$skey\t$count\t$fre1\t$fre\t$hasha{$skey}\t$hashi{$sa}\t$snpid\n";
	}
	print OUT "$key\t$count\t$fre1\t$fre\t$type\t$hash{$key}\t$snpid\n";

}
close OUT;
close COUT;

