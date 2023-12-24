#!/usr/bin/perl -w
my %hash;
open IN, "../extract.near.dbsnp.seq.fa" or die "$!\n";
while(<IN>){
	chomp;
	s/^>//;
	my $id = $_;
	my $seq = <IN>;
	chomp $seq;
	$hash{$id} = $seq;
}

open OUT, ">smSNPs.seq.motif.merge.txt";
print OUT "ID\tmotif\tseq\tCGDHNDHCG\n";
open CI, "sequences.tsv" or die "$!\n";
while(<CI>){
	chomp;
	next if(/^motif/);
	next if(/^#/);
	my ($motif, $id) = (split)[0,3];
	my $seq = "-";
	if(exists $hash{$id}){$seq = $hash{$id}}
	my $tag = (split /-/,$motif)[0];
	if($motif eq "1-NCGDHNDHCGN" && $seq =~ m/(CG\w\w\w\w\wCG)/){
	    $tag = 1.1;
	    my $str = $1;
	    my $str1 = substr($str,2,2);
	    my $str11 = $str1; $str11 =~ tr/ATCG/TAGC/;
	    my $str12 = reverse $str11;
	    my $str13 = reverse $str1;
	    my $str2 = substr($str,5,2);
	    
	    if($str1 eq $str2){$tag=1.2} #CGXYNXYCG
	    
	    if($str11 eq $str2){$tag=1.3} #CGxyNXYCG
	    if($str12 eq $str2){$tag=1.4} #CGyxNXYCG
	    if($str13 eq $str2){$tag=1.5} #CGYXNXYCG
	    #if($id eq "CAMSAP2_1_200833009"){print "$id\t$str\t$str1\t$str2\n"}
    	}
	elsif($motif eq "1-NCGDHNDHCGN" && $seq =~ m/(GC\w\w\w\w\wGC)/){
            $tag = 1.1;
            my $str = $1;
            my $str1 = substr($str,2,2);
            my $str11 = $str1; $str11 =~ tr/ATCG/TAGC/;
            my $str12 = reverse $str11;
            my $str13 = reverse $str1;
            my $str2 = substr($str,5,2);
            
            if($str1 eq $str2){$tag=1.2} #CGXYNXYCG
            
            if($str11 eq $str2){$tag=1.3} #CGxyNXYCG
            if($str12 eq $str2){$tag=1.4} #CGyxNXYCG
            if($str13 eq $str2){$tag=1.5} #CGYXNXYCG
        }

 	print OUT "$id\t$motif\t$seq\t$tag\n";
}
close OUT;

