#!/usr/bin/perl -w

my %hash;
open IN, "extract.near.dbsnp.seq.fa" or die "$!\n";
while(<IN>){
	chomp;
	s/^>//;
	my $id = $_;
	my $seq = <IN>;
	chomp $seq;
	$hash{$id} = $seq; 
}
close IN;

open OUT, ">smSNPs.Palindromes.result.txt";
print OUT "ID\tcontext\tprefix\tsuffix\tseq\tlen\tgap\tmis\n";
foreach my $id(sort keys %hash){
	`../IUPACpal -f extract.near.dbsnp.seq.fa -s $id -m 4 -M 20 -g 5 -x 1 -o tmp/$id.out`;
        my $seq = $hash{$id};
	open TE, "tmp/$id.out" or die "$!\n";
	while(<TE>){
		chomp;
		next if(/^[P|S|M|E|N]/);
		if(/^\d/){
			my ($s1,$nn1,$e1) = (split)[0,1,2];
			my $line2=<TE>; chomp $line2;
			my $line3=<TE>; chomp $line3;
			my $m = $line2 =~ tr/\|/m/;
			my ($s2,$nn2,$e2) = (split /\s+/, $line3)[0,1,2];
			my $str = substr($seq,$s1-1,$s2-$s1+1);
			my $len = length($nn1);
			my $gap = length($str) - 2*$len;
			my $mis = $len - $m;
			print OUT "$id\t$seq\t$nn1\t$nn2\t$str\t$len\t$gap\t$mis\n";
		}
	}
	close TE;
	`rm -rf tmp/$id.out`;
}
close OUT;

