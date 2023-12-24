#!/usr/bin/perl -w

my ($infile, $ofile) = @ARGV;
my (%hash,%hashi,%hash1,%hash2,%hash3,%hashs);

open CI, "../tcga.clin.info.txt" or die "$!\n";
<CI>;
while(<CI>){
        chomp;
	s/\#N\/A/NA/g;
	my ($id,$type,$age,$sex,$stage,$os,$ostime,$pfs,$pfstime) = (split /\t/,$_)[1,2,3,4,6,25,26,31,32];

	$hashi{$id} = "$type\t$age\t$sex\t$stage\t$os\t$ostime\t$pfs\t$pfstime";
}
close CI;

open IN, "select.snp.high.freq.list" or die "$!\n";
while(<IN>){
	chomp;
	my ($gene,$chr,$pos,$cancer) = (split)[0,1,2,3];
	my $k = "chr$chr".":"."$pos";
	$hash{$cancer}{$k} = 1;
}

open OUT1, ">STAD.merge.somatic.snv.maf";
open OUT2, ">UCEC.merge.somatic.snv.maf";
open OUT3, ">COAD.merge.somatic.snv.maf";

open CI, "../TCGA.merge.somatic.snv.maf" or die "$!\n";
<CI>;<CI>;<CI>;<CI>;<CI>;<CI>;<CI>; my $head = <CI>;
print OUT1 "$head";
print OUT2 "$head";
print OUT3 "$head";
while(<CI>){
        chomp;
        next if(/^#/);
        next if(/^Hugo_Symbol/);
        my ($chr,$pos,$ref,$alt,$sa,$tref,$talt,$nref,$nalt) = (split /\t/,$_)[4,5,10,12,15,40,41,43,44];
        my $id = substr($sa,0,12);
	$chr =~ s/chr//;
	next if($nalt >=3);
        my $tsum = $talt + $tref;
        next if($tsum < 5);
        my $tfre = $talt / $tsum;
        next if($tfre < 0.05);

	if(exists $hashi{$id}){
		my $info = $hashi{$id};
		my $cancer = (split /\t/,$info)[0];
		if($cancer eq "STAD"){print OUT1 "$_\n";}
		if($cancer eq "UCEC"){print OUT2 "$_\n";}
		if($cancer eq "COAD"){print OUT3 "$_\n";}	
			
        }
}
close CI;
close OUT1;close OUT2;close OUT3;

my $hh = "Tumor_Sample_Barcode\tSample\tCancer\tAge\tSex\tStage\tOS\tOStime\tPFS\tPFStime";
foreach my $cancer(sort keys %hash){
	my $maffile = "$cancer.merge.somatic.snv.maf";
	my %hashp = %{$hash{$cancer}};
	my %hashm;
	open TI, "$maffile" or die "$!\n";
	<TI>;
	while(<TI>){
		chomp;
		my ($chr,$pos,$sa) = (split)[4,5,15];
		my $ssa = substr($sa,0,12);
		$hashs{$ssa} = $sa;
		my $kk = "$chr".":"."$pos";
		my $key = "$ssa\t$kk";
		if(exists $hashp{$kk}){$hashm{$key} = 1;}
	}
	close TI;
	my @snps;
	foreach my $id(sort keys %hashp){push @snps, $id;}
	open OUT, ">$cancer.clin.info.txt";
	my $Snps = join("\t",@snps);
	print OUT "$hh\t$Snps\n";
	foreach my $id(sort keys %hashi){
		my $ca = (split /\t/,$hashi{$id})[0];
		if($ca eq $cancer){
			my $oo = "";
			foreach my $snp(@snps){
				my $key = "$id\t$snp";
				if(exists $hashm{$key}){$oo .= "\t1"}
				else{$oo .= "\t0";}
			}
			#my $sid = "-"; if(exists $hashs{$id}){$sid = $hashs{$id}}
			next if(!exists $hashs{$id});
			my $sid = $hashs{$id};
			print OUT "$sid\t$id\t$hashi{$id}"."$oo\n";
		}
	}
	close OUT;

}

