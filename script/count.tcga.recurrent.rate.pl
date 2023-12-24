#!/usr/bin/perl
my ($sum1,$count1,$sum2,$count2) = (0,0,0,0);
my (%hash,%hash1,%hash2,%hasha);
open CI, "summary.frequency.txt" or die "$!\n";
while(<CI>){
	chomp;
	my ($gene,$chr,$pos,$count) = (split)[0,1,2,3];
	my $key = "$chr\t$pos";
	$hash{$key} = 1;
	$hash2{$count} += 1;
	$hasha{$count} = 1;
}
close CI;
open OUT, ">tcga.snv.site.no.smSNPs.count.txt";
open IN, "tcga.snv.site.count.txt" or die "$!\n";
while(<IN>){
	chomp;
	s/^\s//g;
	my ($c,$gene,$chr,$pos) = (split)[0,1,2,3];
	$chr =~ s/chr//g;
	my $key = "$chr\t$pos";
	
	if(exists $hash{$key}){
		$hash{$key} += 1;
		$sum2 += 1;
		if($c >1){$count2 += 1}
	}else{
		$sum1 += 1;
		if($c > 1){$count1 +=1}
		print OUT "$c\t$gene\t$chr\t$pos\n";
		$hash1{$c} += 1;
		$hasha{$c} = 1;
	}
}
print "$sum1\t$count1\t$sum2\t$count2\n";
my $fre1 = $count1/$sum1;
my $fre2 = $count2/$sum2;

print "tcga recurrent rate:$fre1\n";
print "smSNPs recurrent rate:$fre2\n";
my $sas = 10574;
open OUT2, ">smSNPs.recurrent.freq.dis.txt";
print OUT2 "Count\tfreq\tnum\tgroup\n";
#print OUT2 "Count\tTCGA\tsmSNPs\tTCGA_num\tsmSNPs_num\n";
foreach my $c(sort {$a <=> $b} keys %hasha){
	my $fre1 = 0;
	my $fre2 = 0;
	if(exists $hash1{$c}){$fre1 = $hash1{$c} / $sum1; print OUT2 "$c\t$fre1\t$hash1{$c}\tTCGA\n"}
	if(exists $hash2{$c}){$fre2 = $hash2{$c} / $sum2; print OUT2 "$c\t$fre2\t$hash2{$c}\tsmSNPs\n"}
	#print OUT2 "$c\t$fre1\t$fre2\t$hash1{$c}\t$hash2{$c}\n";
}

close OUT2;
