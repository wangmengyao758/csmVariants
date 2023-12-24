library(maftools)

vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)


ucec.maf1 <- "UCEC.merge.somatic.snv.chr1:154590148.mut.maf"
ucec.maf2 <- "UCEC.merge.somatic.snv.chr1:154590148.wild.maf"
ucec.clin <- "UCEC.clin.info.txt"

ucec1 = read.maf(maf = ucec.maf1, clinicalData = ucec.clin)
ucec2 = read.maf(maf = ucec.maf2, clinicalData = ucec.clin)
write.mafSummary(maf = ucec1, basename = 'UCEC.ADAR.mut')
write.mafSummary(maf = ucec2, basename = 'UCEC.ADAR.wild')

pdf(file="UCEC.ADAR.compare.pdf")
plotmafSummary(maf = ucec1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = ucec2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = ucec1, top = 30)
oncoplot(maf = ucec2, top = 30)

OncogenicPathways(maf = ucec1)
OncogenicPathways(maf = ucec2)

#ucec1.tnm = trinucleotideMatrix(maf = ucec1, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
#ucec2.tnm = trinucleotideMatrix(maf = ucec2, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

#plotApobecDiff(tnm = ucec1.tnm, maf = ucec1, pVal = 0.2)
#plotApobecDiff(tnm = ucec2.tnm, maf = ucec2, pVal = 0.2)



pt.vs.rt <- mafCompare(m1 = ucec1, m2 = ucec2, m1Name = 'ADAR_mut', m2Name = 'ADAR_wild', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$result, file="UCEC.ADAR.gene.compare.txt",quote=F,sep="\t",col.names=T,row.names=F)

#forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)
#genes <- c("PTEN","ARID1A","PIK3CA","KRAS","PIK3R1","CCND1","CTCF","KMT2B","TP53","FOXD4","RNF43","JAK1","INPPL1","ZFHX3","WBP1")
genes <- c("PTEN","ARID1A","RNF43","PIK3CA","KMT2B","ZFHX3","JAK1","PIK3R1","CTCF","KRAS","INPPL1","WBP1","TP53" )


coOncoplot(m1 = ucec1, m2 = ucec2, m1Name = 'ADAR_mut', m2Name = 'ADAR_wild', genes = genes, removeNonMutated = TRUE)
coBarplot(m1 = ucec1, m2 = ucec2, m1Name = 'ADAR_mut', m2Name = 'ADAR_wild', genes=genes, colors=vc_cols)
dev.off()

ucec.maf1 <- "UCEC.merge.somatic.snv.chr5:88206508.mut.maf"
ucec.maf2 <- "UCEC.merge.somatic.snv.chr5:88206508.wild.maf"
ucec1 = read.maf(maf = ucec.maf1, clinicalData = ucec.clin)
ucec2 = read.maf(maf = ucec.maf2, clinicalData = ucec.clin)
write.mafSummary(maf = ucec1, basename = 'UCEC.TEME161B.mut')
write.mafSummary(maf = ucec2, basename = 'UCEC.TMEM161B.wild')

pdf(file="UCEC.TMEM161B.compare.pdf")
plotmafSummary(maf = ucec1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = ucec2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = ucec1, top = 30)
oncoplot(maf = ucec2, top = 30)
OncogenicPathways(maf = ucec1)
OncogenicPathways(maf = ucec2)

pt.vs.rt <- mafCompare(m1 = ucec1, m2 = ucec2, m1Name = 'TMEM161B_mut', m2Name = 'TMEM161B_wild', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$result, file="UCEC.TMEM161B.gene.compare.txt",quote=F,sep="\t",col.names=T,row.names=F)
#forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)
genes <- c("PTEN","ARID1A","PIK3CA","KRAS","PIK3R1","CCND1","CTCF","UPF3A","KMT2B","DLL1","PRRT2","TP53","FOXD4","RNF43","JAK1","INPPL1","ZFHX3","WBP1","MSH2")
#genes <- c("PTEN","ARID1A","TTN","KMT2D","ZFHX3","KMT2B","JAK1","PIK3CA","MUC16","MUC4","TP53")
coOncoplot(m1 = ucec1, m2 = ucec2, m1Name = 'TMEM161B_mut', m2Name = 'TMEM161B_wild', genes = genes, removeNonMutated = TRUE)
coBarplot(m1 = ucec1, m2 = ucec2, m1Name = 'TMEM161B_mut', m2Name = 'TMEM161B_wild', genes=genes, colors=vc_cols)
dev.off()

