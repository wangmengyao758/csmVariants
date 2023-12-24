library(sigminer)
library(maftools)

stad1 <- "STAD.merge.somatic.snv.chr1:183544857.mut.maf"
stad2 <- "STAD.merge.somatic.snv.chr1:183544857.wild.maf"
stad3 <- "STAD.merge.somatic.snv.chr14:50823156.mut.maf"
stad4 <- "STAD.merge.somatic.snv.chr14:50823156.wild.maf"

ucec1 <- "UCEC.merge.somatic.snv.chr1:154590148.mut.maf"
ucec2 <- "UCEC.merge.somatic.snv.chr1:154590148.wild.maf"
ucec3 <- "UCEC.merge.somatic.snv.chr5:88206508.mut.maf"
ucec4 <- "UCEC.merge.somatic.snv.chr5:88206508.wild.maf"

sdata1 <- read.maf(stad1)
slotNames(sdata1)
mt_tally_all1 <- sig_tally(sdata1,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode="ALL", add_trans_bias = TRUE)
save(mt_tally_all1, file="STAD.chr1:183544857.mut.sigminer.RData")
write.table(mt_tally_all1$APOBEC_scores, file="STAD.chr1:183544857.mut.APOBEC_scores.txt",sep="\t",quote=F,row.names=F)
write.table(mt_tally_all1$SBS_96, file="STAD.chr1:183544857.mut.SBS_96.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all1$ID_83, file="STAD.chr1:183544857.mut.ID_83.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all1$DBS_78, file="STAD.chr1:183544857.mut.DBS_78.txt",sep="\t",quote=F,row.names=T)
pdf(file="show_catalogue.STAD.chr1.183544857.mut.pdf",height=4,width=10)
show_catalogue(t(mt_tally_all1$SBS_96), mode="SBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all1$SBS_96))
show_catalogue(t(mt_tally_all1$ID_83), mode="ID", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all1$ID_83))
show_catalogue(t(mt_tally_all1$DBS_78), mode="DBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all1$DBS_78))
dev.off()

sbs_sig_fit <- sig_fit(t(mt_tally_all1$SBS_96),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit, file="STAD.chr1:183544857.mut.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit <- sig_fit(t(mt_tally_all1$ID_83),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit, file="STAD.chr1:183544857.mut.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit <- sig_fit(t(mt_tally_all1$DBS_78),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit, file="STAD.chr1:183544857.mut.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)

sdata2 <- read.maf(stad2)
slotNames(sdata2)
mt_tally_all2 <- sig_tally(sdata2,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode="ALL", add_trans_bias = TRUE)
save(mt_tally_all2, file="STAD.chr1:183544857.wild.sigminer.RData")
write.table(mt_tally_all2$APOBEC_scores, file="STAD.chr1:183544857.wild.APOBEC_scores.txt",sep="\t",quote=F,row.names=F)
write.table(mt_tally_all2$SBS_96, file="STAD.chr1:183544857.wild.SBS_96.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all2$ID_83, file="STAD.chr1:183544857.wild.ID_83.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all2$DBS_78, file="STAD.chr1:183544857.wild.DBS_78.txt",sep="\t",quote=F,row.names=T)
pdf(file="show_catalogue.STAD.chr1.183544857.wild.pdf",height=4,width=10)
show_catalogue(t(mt_tally_all2$SBS_96), mode="SBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all2$SBS_96))
show_catalogue(t(mt_tally_all2$ID_83), mode="ID", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all2$ID_83))
show_catalogue(t(mt_tally_all2$DBS_78), mode="DBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all2$DBS_78))
dev.off()

sbs_sig_fit <- sig_fit(t(mt_tally_all2$SBS_96),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit, file="STAD.chr1:183544857.wild.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit <- sig_fit(t(mt_tally_all2$ID_83),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit, file="STAD.chr1:183544857.wild.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit <- sig_fit(t(mt_tally_all2$DBS_78),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit, file="STAD.chr1:183544857.wild.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)


sdata3 <- read.maf(stad3)
slotNames(sdata3)
mt_tally_all3 <- sig_tally(sdata3,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode="ALL", add_trans_bias = TRUE)
save(mt_tally_all3, file="STAD.chr14:50823156.mut.sigminer.RData")
write.table(mt_tally_all3$APOBEC_scores, file="STAD.chr14:50823156.mut.APOBEC_scores.txt",sep="\t",quote=F,row.names=F)
write.table(mt_tally_all3$SBS_96, file="STAD.chr14:50823156.mut.SBS_96.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all3$ID_83, file="STAD.chr14:50823156.mut.ID_83.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all3$DBS_78, file="STAD.chr14:50823156.mut.DBS_78.txt",sep="\t",quote=F,row.names=T)
pdf(file="show_catalogue.STAD.chr14.50823156.mut.pdf",height=4,width=10)
show_catalogue(t(mt_tally_all3$SBS_96), mode="SBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all3$SBS_96))
show_catalogue(t(mt_tally_all3$ID_83), mode="ID", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all3$ID_83))
show_catalogue(t(mt_tally_all3$DBS_78), mode="DBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all3$DBS_78))
dev.off()

sbs_sig_fit <- sig_fit(t(mt_tally_all3$SBS_96),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit, file="STAD.chr14:50823156.mut.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit <- sig_fit(t(mt_tally_all3$ID_83),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit, file="STAD.chr14:50823156.mut.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit <- sig_fit(t(mt_tally_all3$DBS_78),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit, file="STAD.chr14:50823156.mut.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)


sdata4 <- read.maf(stad4)
slotNames(sdata4)
mt_tally_all4 <- sig_tally(sdata4,ref_genome = "BSgenome.Hsapiens.UCSC.hg38",use_syn = TRUE,mode="ALL", add_trans_bias = TRUE)
save(mt_tally_all4, file="STAD.chr14:50823156.wild.sigminer.RData")
write.table(mt_tally_all4$APOBEC_scores, file="STAD.chr14:50823156.wild.APOBEC_scores.txt",sep="\t",quote=F,row.names=F)
write.table(mt_tally_all4$SBS_96, file="STAD.chr14:50823156.wild.SBS_96.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all4$ID_83, file="STAD.chr14:50823156.wild.ID_83.txt",sep="\t",quote=F,row.names=T)
write.table(mt_tally_all4$DBS_78, file="STAD.chr14:50823156.wild.DBS_78.txt",sep="\t",quote=F,row.names=T)
pdf(file="show_catalogue.STAD.chr14.50823156.wild.pdf",height=4,width=10)
show_catalogue(t(mt_tally_all4$SBS_96), mode="SBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all4$SBS_96))
show_catalogue(t(mt_tally_all4$ID_83), mode="ID", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all4$ID_83))
show_catalogue(t(mt_tally_all4$DBS_78), mode="DBS", style = "cosmic", x_label_angle = 90, y_tr = function(x) x / sum(mt_tally_all4$DBS_78))
dev.off()

sbs_sig_fit <- sig_fit(t(mt_tally_all4$SBS_96),sig_index="ALL", sig_db="SBS", mode="SBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(sbs_sig_fit, file="STAD.chr14:50823156.wild.sbs_sig_fit.txt",sep="\t",quote=F,row.names=F)
id_sig_fit <- sig_fit(t(mt_tally_all4$ID_83),sig_index="ALL", sig_db="ID", mode="ID", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(id_sig_fit, file="STAD.chr14:50823156.wild.id_sig_fit.txt",sep="\t",quote=F,row.names=F)
dbs_sig_fit <- sig_fit(t(mt_tally_all4$DBS_78),sig_index="ALL", sig_db="DBS", mode="DBS", return_class = "data.table", rel_threshold = 0.05,type="relative")
write.table(dbs_sig_fit, file="STAD.chr14:50823156.wild.dbs_sig_fit.txt",sep="\t",quote=F,row.names=F)

