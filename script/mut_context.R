args=commandArgs(T)
vcffile<-args[2]
sample<-args[1]

outfile <- paste(sample,".context.txt",sep="")
matfile <- paste(sample,".96context.txt",sep="")
pdf1 <- paste(sample, ".spectrum.pdf",sep="")
pdf2 <- paste(sample, ".96.profile.pdf",sep="")

library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

vcfs <- read_vcfs_as_granges(vcffile, sample, ref_genome)

context <- mut_context(vcfs[[1]], ref_genome)
sdata <- as.data.frame(vcfs[[1]])
sdata$context <- context
sdata$ALT <- NULL
#sdata$REF <- NULL

write.table(sdata,outfile,col.names=T,row.names=T,sep="\t",quote=F)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
print(type_occurrences)

pdf(pdf1,height=5,width=10)
plot_spectrum(type_occurrences, CT = TRUE,error_bars = 'none')
dev.off()

#mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
write.table(mut_mat,matfile,col.names=T,row.names=T,sep="\t",quote=F)
pdf(pdf2,height=3,width=10)
plot_96_profile(mut_mat)
dev.off()
