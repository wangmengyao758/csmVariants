#!/usr/bin/perl -w

open IN, "tcga.Palindromes.result.txt" or die "$!\n";
open OUT, ">tcga.Palindromes.result.filter.txt";
my $head = <IN>; chomp $head;
print OUT "$head\tstart\tregion\n";
my (%hash, %hashs, %hashi);
while(<IN>){
	chomp;
	my ($id,$fa,$seq,$len,$gap,$mis) = (split)[0,1,4,5,6,7];
        my $key = "$id\t$seq";
	$hash{$key} = $_;
	next if($len<=4 && $mis==1);
	my $score = $len - $mis/2;
	$fa =~ m/$seq/g;
	my $l = length($seq);
	my $e = pos $fa ;
	my $s = $e - $l + 1;
	my $tag = "-"; # mutation on other region;
        my $i2 = $s + $len;
	my $i3 = $s + $len + $gap;
	my $i4 = $s + $l;
	if($s<=11 && $i2>=11){$tag = "flank"}
	if($i2<11 && $i3>11){$tag = "loop"}
	if($i3<=11 && $i4>=11){$tag = "flank"}


	if(!exists $hashs{$id}){$hashs{$id} = $score; $hashi{$id} = "$_\t$s\t$tag"}
	elsif($score > $hashs{$id}){$hashs{$id} = $score; $hashi{$id} = "$_\t$s\t$tag";}

}
close IN;

foreach my $key(sort keys %hashi){
	print OUT "$hashi{$key}\n";
}
close OUT;

