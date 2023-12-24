library(sigminer)
library(maftools)

maffile <- "../TCGA.merge.somatic.snv.dbsnp.maf"
sdata <- read.maf(maffile)
slotNames(sdata)
mt_tally <- sig_tally(sdata,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE)

mt_tally_all <- sig_tally(sdata,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode="ALL", add_trans_bias = TRUE)
#mt_tally_DBS <- sig_tally(sdata,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode = "DBS",add_trans_bias = TRUE)
#mt_tally_ID <- sig_tally(sdata,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode = "ID",add_trans_bias = TRUE)
ls(mt_tally_all)
save(mt_tally_all, file="smSNPs.sigminer.RData")
write.table(mt_tally_all$APOBEC_scores, file="smSNPs.APOBEC_scores.txt",sep="\t",quote=F,row.names=F)
write.table(mt_tally_all$SBS_96, file="smSNPs.SBS_96.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all$ID_83, file="smSNPs.ID_83.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all$DBS_78, file="smSNPs.DBS_78.txt",sep="\t",quote=F,row.names=T)
pdf(file="show_catalogue.smSNPs.pdf",height=4,width=10)
show_catalogue(t(mt_tally_all$SBS_96), mode="SBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all$SBS_96))
show_catalogue(t(mt_tally_all$ID_83), mode="ID", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all$ID_83))
show_catalogue(t(mt_tally_all$DBS_78), mode="DBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all$DBS_78))
dev.off()

sbs_sig_fit <- sig_fit(t(mt_tally_all$SBS_96),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit, file="smSNPs.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit <- sig_fit(t(mt_tally_all$ID_83),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit, file="smSNPs.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit <- sig_fit(t(mt_tally_all$DBS_78),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit, file="smSNPs.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)

#perl merge.cancertype.pl
sbs <- read.table("smSNPs.cancertype.SBS_96.txt")
id <- read.table("smSNPs.cancertype.ID_83.txt")
dbs <- read.table("smSNPs.cancertype.DBS_78.txt")
colnames(sbs) <- colnames(mt_tally_all$SBS_96)
colnames(id) <- colnames(mt_tally_all$ID_83)
colnames(dbs) <- colnames(mt_tally_all$DBS_78)

sbs_sig_fit2 <- sig_fit(t(sbs),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit2, file="smSNPs.cancertype.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit2 <- sig_fit(t(id),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit2, file="smSNPs.cancertype.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit2 <- sig_fit(t(dbs),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit2, file="smSNPs.cancertype.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)








#mt_sig2 <- sig_auto_extract(mt_tally_all$nmf_matrix, K0 = 10, nrun = 10, strategy = "stable")
#sim <- get_sig_similarity(mt_sig2)
#pheatmap::pheatmap(sim$similarity)
#sim_v3 <- get_sig_similarity(mt_sig2, sig_db = "SBS")
#i <- which.max(apply(mt_tally$nmf_matrix, 1, sum))
#example_mat <- mt_tally$nmf_matrix[i, , drop = FALSE] %>% t()
#sig_fit(example_mat, sig_index = 1:30,type = "relative")
#sig_fit(t(mt_tally$nmf_matrix[1:5, ]), sig_index = 1:30, return_class = "data.table", rel_threshold = 0.05,type="relative")

#sigprofiler_extract(mt_tally_all$nmf_matrix, "SigProfiler_smSNPs",  use_conda = FALSE, py_path = "/home/wangmengyao/anaconda3/envs/py36/bin/python")



