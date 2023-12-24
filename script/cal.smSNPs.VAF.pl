#!/usr/bin/perl -w

my (%hash, %hash1,%hash2);
open CI, "../summary.frequency.txt" or die "$!\n";
while(<CI>){
        chomp;
        my ($gene,$chr,$pos,$count) = (split)[0,1,2,3];
        my $key = "$chr\t$pos";
        $hash{$key} = 1;
	#$hash2{$count} += 1;
	#$hasha{$count} = 1;
}
close CI;

open OUT, ">TCGA.readcount.txt";
print OUT "Class\tChr\tPos\tRef\tAlt\tType\tSample\ttref\ttalt\tnref\tnalt\ttvaf\tnvaf\n";
open CI, "../TCGA.merge.somatic.snv.maf" or die "$!\n";
while(<CI>){
        chomp;
        next if(/^#/);
        next if(/^Hugo_Symbol/);
        my ($chr,$pos,$type,$ref,$alt,$sa,$tref,$talt,$nsum,$nref,$nalt) = (split /\t/,$_)[4,5,9,10,12,15,40,41,42,43,44];
        my $id = substr($sa,0,12);
        $chr =~ s/chr//;
        my $key = "$chr\t$pos";
	next if($nalt >=3);
        my $tsum = $talt + $tref;
        next if($tsum < 5);
	next if($nsum < 3);
        my $tfre = $talt / $tsum;
        my $nfre = $nalt / $nsum;
	next if($tfre < 0.05);
	my $class = "TCGA";
	if(exists $hash{$key}){$class = "smSNPs"}
        print OUT "$class\t$chr\t$pos\t$ref\t$alt\t$type\t$sa\t$tref\t$talt\t$nref\t$nalt\t$tfre\t$nfre\n";

}
close CI;
close OUT;

