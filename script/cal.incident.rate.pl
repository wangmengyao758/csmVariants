#!/usr/bin/perl

my $common_loci = 0; #common snp loci in exon region
my $non_common_loci = 0; #non-common snp loci in exon region total_exon - common_loci
my $other_mutation = 0 ;  #exonic non-smSNPs somatic mutation count;
my $total_exon = 0; # total base count in exon region
my $smSNP_exon = 0; #exonic smSNP in TCGA
my (%hash,%hashc,%hasha,$hashs);
#TCGA accessiable region "exon"

open IN, "Agilent_Human_Exon_v5_UTRs_Regions" or die "$!\n";
while(<IN>){
    chomp;
    my ($chr,$s,$e) = (split);
    my $region = "$s\t$e";
    my $len = $e - $s + 1;
    $total_exon += $len;
    push(@{$hash{$chr}}, $region);
}
close IN;
print "$total_exon\n";
sub withinExome{
    my ($pos,@regions) = @_;
    foreach my $region(@regions){
        my ($start,$end) = (split /\t/,$region)[0,1];
        if($pos >= $start && $pos <= $end){
            return 1;
        }
    }
    return 0;
}

open IN1, "/home/wangmengyao/data3/dbsnp/hg38/common_all_20180418.vcf" or die "$!\n";
while(<IN1>){
    chomp;
    next if(/^#/);
    my ($chr,$pos) = (split)[0,1];
    $chr = "chr"."$chr";
    my $re = "$chr\t$pos";
    $hash{$re} += 1;
    next if($hash{$re} > 1);
    next if(!exists $hash{$chr});
    my @exon_intervals = @{$hash{$chr}};
    if(withinExome($pos,@exon_intervals)){$common_loci += 1;$hasha{$re}=1}
}
close IN1;
print "$common_loci\n";
open IN2, "/home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.dbsnp.maf" or die "$!\n";
<IN2>;
while(<IN2>){
    chomp;
    my ($chr,$pos) = (split)[4,5];
    my $re = "$chr\t$pos";
    if(exists $hasha{$re}){$smSNPs_exon += 1}
    $hashs{$re} = 1;
}
close IN2;

open IN3, "/home/wangmengyao/data3/dbsnp/TCGA.merge.somatic.snv.maf" or die "$!\n";
while(<IN3>){
    chomp;
    next if(/^#|Hugo_Symbol/);
    my ($chr,$pos) = (split)[4,5];
    my $re = "$chr\t$pos";
    next if(exists $hashs{$re}); ##smSNPs;
    next if(exists $hasht{$re});
    my @exon_intervals = @{$hash{$chr}};
    if(withinExome($pos,@exon_intervals)){$other_mutation += 1;}
}
close IN3;
my $non_common_loci = $total_exon - $common_loci;
print "common snp loci in exon region: $common_loci\n";
print "non-common snp loci in exon region: $non_common_loci\n";
print "exonic smSNP count in TCGA: $smSNP_exon\n";
print "exonic non-smSNPs somatic mutation count: $other_mutation\n";

