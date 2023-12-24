less TCGA.merge.somatic.snv.maf.gz | cut -f 1,5,6 | sort | uniq -c > tcga.snv.site.count.txt

less TCGA.merge.somatic.snv.maf |grep -v "#"|grep -v Hugo_Symbol|cut -f 14|sort|uniq -c >snp.count.txt
less snp.count.txt|grep rs|wc -l 2418383 : snp count
less snp.count.txt|grep rs|awk '$1>10'|wc -l 56114 : high recurrent snp count
less snp.count.txt|grep rs|awk '{sum += $1} END{print sum}' 5172369 : total snp number
less TCGA.merge.somatic.snv.maf |grep -v "#"|grep -v Hugo_Symbol|wc -l  13774739 : total mutation number


perl cal.incident.rate.pl

#CMplot
perl class.pl
Rscript fisher_test.R
perl add.snp.pvalue.pl
perl add.mut.type.pl

#mutational spectrum
perl get.lego.plot.input.pl
Rscript mut_context.R tcga tcga.context.hg38.vcf
Rscript mut_context.R smsnps smsnps.context.hg38.vcf
Rscript mut_context.R dbsnps dbsnps.context.hg38.vcf
library(sigminer)
library(MutationalPatterns)
x1 <- get_sig_db("SBS_hg38")
aa <- x1$db
SBS <- aa[order(row.names(aa)),]

sdata <- read.table("merge.96.context2.txt",header=T,row.names=1)
#smSNPs vs. SBS1
cos_sim(SBS[,1],sdata[order(row.names(sdata)),1])

#mutational signature comparison

Rscript run.sigminer.tcga.R
Rscript run.sigminer.R
Rscript plot.sig.R


#motif enrichment
perl extract.near.seq.pl
/home/wangmengyao/software/meme-5.5.0/bin/streme --p extract.near.dbsnp.seq.fa --o streme_out --dna 
/home/wangmengyao/software/meme-5.5.0/bin/seq --p extract.near.dbsnp.seq.fa --m meme_out1/meme.txt --o sea_out
perl merge.motif.seq.pl

perl run.smSNP.IUPACpal.pl
perl run.tcga.IUPACpal.pl
perl filter.smSNP.IUPACpal.pl
perl filter.tcga.IUPACpal.pl
perl annotation.mutation.Palindromes.pl

#Molecular subtyping
perl subset.maf.pl
perl split.maf.pl STAD chr14:50823156\chr1:183544857\chr2:91699879\chr3:149968578\chr3:27720073\chr8:112503980\chr9:22451040
perl split.maf.pl UCEC chr10:116636766/chr14:50823156/chr1:154590148/chr1:183544857/chr5:88206508/chr6:18258069
Rscript read.maf.R  #STAD
Rscript read.maf2.R   #UCEC
Rscript run.sigminer.subset.R



