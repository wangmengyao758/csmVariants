#!/usr/bin/perl -w
my (%hashi,%hashc,%hash1,%hash2,%hash3,%hash);
open CI, "../tcga.clin.info.txt" or die "$!\n";
<CI>;
while(<CI>){
        chomp;
        my ($id,$cancer,$age,$sex) = (split)[1,2,3,4];
        $hash{$id} = $cancer;
	$hashc{$cancer} += 1;
}
close CI;

open IN1, "tcga.SBS_96.txt" or die "$!\n";
open OUT1, ">tcga.cancertype.SBS_96.txt";
my $header = <IN1>;
chomp $header;
print OUT1 "$header\n";
my @types = split /\t/, $header;
while(<IN1>){
	chomp;
	my @temps = split;
	my $sa = shift @temps;
	my $id = substr($sa,0,12);
	for(my $i=0;$i<=$#types;$i++){
		my $type = $types[$i];
		my $count = $temps[$i];
		next if(!exists $hash{$id});
		my $cancer = $hash{$id};
		my $key = "$cancer\t$type";
		$hash1{$key} += $count;
	}
}
close IN1;
foreach my $c(sort keys %hashc){
	my $out = $c;
	foreach my $type(@types){
		my $key = "$c\t$type";
		my $count = $hash1{$key};
		$out .= "\t$count";
	}
	print OUT1 "$out\n";
}
close OUT1;

open IN2, "tcga.DBS_78.txt" or die "$!\n";
open OUT2, ">tcga.cancertype.DBS_78.txt";
my $header2 = <IN2>;
chomp $header2;
print OUT2 "$header2\n";
my @types2 = split /\t/, $header2;
while(<IN2>){
        chomp;
        my @temps = split;
        my $sa = shift @temps;
        my $id = substr($sa,0,12);
        for(my $i=0;$i<=$#types2;$i++){
                my $type = $types2[$i];
                my $count = $temps[$i];
                next if(!exists $hash{$id});
                my $cancer = $hash{$id};
                my $key = "$cancer\t$type";
                $hash2{$key} += $count;
        }
}
close IN2;
foreach my $c(sort keys %hashc){
        my $out = $c;
        foreach my $type(@types2){
                my $key = "$c\t$type";
                my $count = $hash2{$key};
                $out .= "\t$count";
        }
        print OUT2 "$out\n";
}
close OUT2;

open IN3, "tcga.ID_83.txt" or die "$!\n";
open OUT3, ">tcga.cancertype.ID_83.txt";
my $header3 = <IN3>;
chomp $header3;
print OUT3 "$header3\n";
my @types3 = split /\t/, $header3;
while(<IN3>){
        chomp;
        my @temps = split;
        my $sa = shift @temps;
        my $id = substr($sa,0,12);
        for(my $i=0;$i<=$#types3;$i++){
                my $type = $types3[$i];
                my $count = $temps[$i];
                next if(!exists $hash{$id});
                my $cancer = $hash{$id};
                my $key = "$cancer\t$type";
                $hash3{$key} += $count;
        }
}
close IN3;
foreach my $c(sort keys %hashc){
        my $out = $c;
        foreach my $type(@types3){
                my $key = "$c\t$type";
                my $count = $hash3{$key};
                $out .= "\t$count";
        }
        print OUT3 "$out\n";
}
close OUT3;
