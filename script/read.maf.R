library(maftools)

stad.maf1 <- "STAD.merge.somatic.snv.chr1:183544857.mut.maf"
stad.maf2 <- "STAD.merge.somatic.snv.chr1:183544857.wild.maf"
stad.clin <- "STAD.clin.info.txt"

stad1 = read.maf(maf = stad.maf1, clinicalData = stad.clin)
stad2 = read.maf(maf = stad.maf2, clinicalData = stad.clin)
write.mafSummary(maf = stad1, basename = 'STAD.SMG7.mut')
write.mafSummary(maf = stad2, basename = 'STAD.SMG7.wild')

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

pdf(file="STAD.SMG7.compare.pdf")
plotmafSummary(maf = stad1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = stad2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = stad1, top = 30)
oncoplot(maf = stad2, top = 30)

OncogenicPathways(maf = stad1)
OncogenicPathways(maf = stad2)

#stad1.tnm = trinucleotideMatrix(maf = stad1, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
#stad2.tnm = trinucleotideMatrix(maf = stad2, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

#plotApobecDiff(tnm = stad1.tnm, maf = stad1, pVal = 0.2)
#plotApobecDiff(tnm = stad2.tnm, maf = stad2, pVal = 0.2)



pt.vs.rt <- mafCompare(m1 = stad1, m2 = stad2, m1Name = 'SMG7_mut', m2Name = 'SMG7_wild', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$result, file="STAD.SMG7.gene.compare.txt",quote=F,sep="\t",col.names=T,row.names=F)

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)

#genes = c("TTN", "TP53", "MUC16","ARID1A","OBSCN","KMT2D","KMT2C","MSH3", "MSH6","PIK3CA")
#genes = c("TTN","ACVR2A", "ARID1A", "RPL22", "TTK", "OBSCN","KMT2D","RNF43","UBR5","SYNE1","MSH3","KMT2C","PIK3CA","TP53","MSH6")
genes = c("ARID1A","RPL22","TP53","ERBB3","FBXW7","PIK3CA","BCOR","KRAS","PTEN","B2M","MAP2K7")
coOncoplot(m1 = stad1, m2 = stad2, m1Name = 'SMG7_mut', m2Name = 'SMG7_wild', genes = genes, removeNonMutated = TRUE)
coBarplot(m1 = stad1, m2 = stad2, m1Name = 'SMG7_mut', m2Name = 'SMG7_wild', genes=genes, colors=vc_cols)
dev.off()



stad.maf1 <- "STAD.merge.somatic.snv.chr14:50823156.mut.maf"
stad.maf2 <- "STAD.merge.somatic.snv.chr14:50823156.wild.maf"

stad1 = read.maf(maf = stad.maf1, clinicalData = stad.clin)
stad2 = read.maf(maf = stad.maf2, clinicalData = stad.clin)
write.mafSummary(maf = stad1, basename = 'STAD.NIN.mut')
write.mafSummary(maf = stad2, basename = 'STAD.NIN.wild')
pdf(file="STAD.NIN.compare.pdf")
plotmafSummary(maf = stad1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = stad2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = stad1, top = 30)
oncoplot(maf = stad2, top = 30)

OncogenicPathways(maf = stad1)
OncogenicPathways(maf = stad2)

pt.vs.rt <- mafCompare(m1 = stad1, m2 = stad2, m1Name = 'NIN_mut', m2Name = 'NIN_wild', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$result, file="STAD.NIN.gene.compare.txt",quote=F,sep="\t",col.names=T,row.names=F)

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.1)


#genes = c("TTN", "TP53", "MUC16","ARID1A","PLEC","OBSCN","KMT2D","KMT2C","MSH3", "MSH6","PIK3CA")
#genes = c("TTN", "ARID1A", "KMT2D","OBSCN","PLEC","MUC16","KMT2C","MSH3","TP53","MSH6","PIK3CA")
coOncoplot(m1 = stad1, m2 = stad2, m1Name = 'NIN_mut', m2Name = 'NIN_wild', genes = genes, removeNonMutated = TRUE)
coBarplot(m1 = stad1, m2 = stad2, m1Name = 'NIN_mut', m2Name = 'NIN_wild', genes=genes, colors=vc_cols)
dev.off()

