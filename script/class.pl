#!/usr/bin/perl -w

my @cancers = ("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM");
my @sexs = ("MALE", "FEMALE");
my @ages = ("<60",">=60");
my @tissues;
my @types;

sub uniq {my %seen; return grep {!$seen{$_}++} @_;}

my (%hashe, %hash,%hashc,%hashs,%hasha,%hasht,%hashb);
open TI, "tcga.clin.info.txt" or die "$!\n";
<TI>;
while(<TI>){
	chomp;
	my $type = (split)[2];
	$hashe{$type} +=1;
}
close TI;
open IN, "summary.frequency.table.txt" or die "$!\n";
while(<IN>){
	chomp;
	my ($sa,$gene,$chr,$pos,$count,$fret,$freq,$type,$cancer,$age,$gender) = (split)[0,1,2,3,4,5,6,7,8,9,10];
	my $key = "$gene\t$chr\t$pos";
        $hash{$key} = "$count\t$fret\t$freq";
	my $ssa = (split /-/, $sa)[3];
	$ssa =~ s/\w$//;
	next if(!$cancer);
	my $k1 = "$key\t$cancer";
	my $Age = "<60";
	if($age =~ /\#/){$Age = "NA"}
	elsif($age >= 60){$Age = ">=60"}
	my $k2 = "$key\t$Age";
	my $k3 = "$key\t$gender";
        $hashc{$k1} += 1;
	$hashs{$k3} += 1;
	$hasha{$k2} += 1;
	my $k4 = "$key\t$type";
	my $k5 = "$key\t$ssa";
	$hasht{$k4} += 1;
	$hashb{$k5} += 1;
	push @tissues, $ssa;
	push @types, $type;
}
close IN;
 @tissues = uniq(@tissues);
 @types = uniq(@types);

 open OUT, ">summary.tcga.cancertype.dbsnp.table.txt";
 my $head1 = "Gene\tChr\tPos\tCount\tFrequency\tdbSNP_freq\tMale\tFemale\tage1\tage2\ttype";
 my $head2 = join("\t",@cancers);
 my $head3 = join("\t",@types);
 my $head4 = join("\t",@tissues);
 $head4 =~ s/0/T0/g;
 $head3 =~ s/\'/_/g;
 print OUT "$head1\t$head2\t$head3\t$head4\n";
 foreach my $key(sort keys %hash){
      my $out1 = $hash{$key};
      my $count = (split /\t/,$out1)[0];
      next if($count < 10);
      my ($out2, $out3, $out4);
      my ($male, $female, $age1, $age2) = (0,0,0,0);
      if(exists $hashs{"$key\tMALE"}){$male = sprintf "%.5f", $hashs{"$key\tMALE"} / $count }
      if(exists $hashs{"$key\tFEMALE"}){$female = sprintf "%.5f", $hashs{"$key\tFEMALE"} / $count }
      if(exists $hasha{"$key\t<60"}){$age1 = sprintf "%.5f", $hasha{"$key\t<60"} / $count }
      if(exists $hasha{"$key\t>=60"}){$age2 = sprintf "%.5f", $hasha{"$key\t>=60"}/ $count}
      foreach my $cancer(@cancers){
	     my $tt = 0;
	     my $k1 = "$key\t$cancer";
	     my $sum = $hashe{$cancer};
	     if(exists $hashc{$k1}){$tt = sprintf "%.5f", $hashc{$k1}/ $sum}
 	     $out2 .= "\t$tt";
      } 
      my ($mut_type,$mscore) = ("",0);
      foreach my $type(@types){
	      my $tt = 0;
	      my $k = "$key\t$type";
	      if(exists $hasht{$k}){$tt = sprintf "%.5f", $hasht{$k} / $count}
	      if($tt>$mscore){$mscore=$tt;$mut_type=$type}
	      $out3 .= "\t$tt";
      }
      foreach my $tiss(@tissues){
	      my $tt = 0;
	      my $k = "$key\t$tiss";
	      if(exists $hashb{$k}){$tt = sprintf "%.5f", $hashb{$k} / $count}
	      $out4 .= "\t$tt";
      }
      print OUT "$key\t$out1\t$male\t$female\t$age1\t$age2\t$mut_type"."$out2"."$out3"."$out4\n"; 
}
 close OUT;


