library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
colors <- rev(pal_npg(alpha =0.9)(2))

sdata <- read.table("smSNPs.cancertype.sbs_sig_fit.txt",header=T,sep="\t",row.names=1)
sdata1 <- sdata[colSums(sdata)>1]
sdata1$SBS54 <- NULL
idata <- read.table("smSNPs.cancertype.id_sig_fit.txt",header=T,sep="\t",row.names=1)
idata1 <- idata[colSums(idata)>1]

ddata <- read.table("smSNPs.cancertype.dbs_sig_fit.txt",header=T,sep="\t",row.names=1)
ddata1 <- ddata[colSums(ddata)>1]

idata1$sample <- rownames(idata1)
idata2 <- melt(idata1)
idata2$group <- "smSNPs"

sdata1$sample <- rownames(sdata1)
sdata2 <- melt(sdata1)
sdata2$group <- "smSNPs"

ddata1$sample <- rownames(ddata1)
ddata2 <- melt(ddata1)
ddata2$group <- "smSNPs"

stdata <- read.table("tcga.cancertype.sbs_sig_fit.txt",header=T,sep="\t",row.names=1)
stdata$sample <- rownames(stdata)
stdata1 <- stdata[,colnames(sdata1)]

itdata <- read.table("tcga.cancertype.id_sig_fit.txt",header=T,sep="\t",row.names=1)
itdata$sample <- rownames(itdata)
itdata1 <- itdata[,colnames(idata1)]

dtdata <- read.table("tcga.cancertype.dbs_sig_fit.txt",header=T,sep="\t",row.names=1)
dtdata$sample <- rownames(dtdata)
dtdata1 <- dtdata[,colnames(ddata1)]

stdata2 <- melt(stdata1)
stdata2$group <- "TCGA"
itdata2 <- melt(itdata1)
itdata2$group <- "TCGA"
dtdata2 <- melt(dtdata1)
dtdata2$group <- "TCGA"

msdata <- rbind(sdata2, stdata2)
midata <- rbind(idata2, itdata2)
mddata <- rbind(ddata2, dtdata2)
colnames(msdata) <- c("cancer","signature","contribution","group")
colnames(midata) <- c("cancer","signature","contribution","group")
colnames(mddata) <- c("cancer","signature","contribution","group")


pdf(file="SBS.signature.compare.cancertype.pdf",height=6,width=10)
ggplot(data = msdata, aes(x=cancer,y=contribution,fill=cancer)) + geom_bar(stat="identity",show.legend = F) + facet_grid(signature ~ group) +scale_y_continuous() + theme_bw() + theme(axis.text.x=element_text(size=8,angle =90,hjust =0.5,vjust = 0.5), axis.text.y = element_text(size=8), strip.text.y = element_text(size=10), panel.grid = element_blank()) + xlab("Cancer Type") + ylab("contribution")
dev.off()

pdf(file="DBS.signature.compare.cancertype.pdf",height=6,width=10)
ggplot(data = mddata, aes(x=cancer,y=contribution,fill=cancer)) + geom_bar(stat="identity",show.legend = F) + facet_grid(signature ~ group) +scale_y_continuous() + theme_bw() + theme(axis.text.x=element_text(size=8,angle =90,hjust =0.5,vjust = 0.5), axis.text.y = element_text(size=8), strip.text.y = element_text(size=10), panel.grid = element_blank()) + xlab("Cancer Type") + ylab("contribution")
dev.off()

pdf(file="ID.signature.compare.cancertype.pdf",height=6,width=10)
ggplot(data = midata, aes(x=cancer,y=contribution,fill=cancer)) + geom_bar(stat="identity",show.legend = F) + facet_grid(signature ~ group) +scale_y_continuous() + theme_bw() + theme(axis.text.x=element_text(size=8,angle =90,hjust =0.5,vjust = 0.5), axis.text.y = element_text(size=8), strip.text.y = element_text(size=10), panel.grid = element_blank()) + xlab("Cancer Type") + ylab("contribution")
dev.off()

sigs <- c("SBS1", "SBS10a", "SBS10b", "SBS15")
sdata <- read.table("merge.SBS.signature.cosmic.contribution.txt",header=T,sep="\t")
ssdata <- sdata[sdata$smSNPs_num>500 & sdata$Cancer == "UCEC",]
ssdata1 <- ssdata[ssdata$signature %in% sigs,]

p1 <- ggplot(ssdata1,aes(x=signature,y=Sample,fill=Sig_smSNPs)) + xlab("") + ylab("") + geom_tile(color="white",size=0.1) + scale_fill_gradient(low = "gray95", high = "tomato") + geom_text(aes(label=round(Sig_smSNPs,2)),angle=45, size=2) +  ggtitle("smSNPs") + theme_minimal() +  theme(panel.grid = element_blank(), axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5, size=8))

p2 <- ggplot(ssdata1,aes(x=signature,y=Sample,fill=Sig_TCGA)) + xlab("") + ylab("") + geom_tile(color="white",size=0.1) + scale_fill_gradient(low = "gray95", high = "tomato") + geom_text(aes(label=round(Sig_TCGA,2)),angle=45, size=2) +  ggtitle("TCGA") + theme_minimal() +  theme(panel.grid = element_blank(), axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5, size=8))

pp <- ggarrange(p1,p2,ncol=2,nrow=1,common.legend=T,legend="bottom",labels=c("A","B"))
pdf(file="merge.SBS.signature.smSNPs.compare.pdf", height=8, width=10)
pp
dev.off()


sigs <- c("ID2", "ID9", "ID14", "ID12")
sdata <- read.table("merge.ID.signature.cosmic.contribution.txt",header=T,sep="\t")
ssdata <- sdata[sdata$smSNPs_num>300 & sdata$Cancer == "UCEC",]
ssdata1 <- ssdata[ssdata$signature %in% sigs,]

p1 <- ggplot(ssdata1,aes(x=signature,y=Sample,fill=Sig_smSNPs)) + xlab("") + ylab("") + geom_tile(color="white",size=0.1) + scale_fill_gradient(low = "gray95", high = "tomato") + geom_text(aes(label=round(Sig_smSNPs,2)),angle=45, size=2) +  ggtitle("smSNPs") + theme_minimal() +  theme(panel.grid = element_blank(), axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5, size=8))

p2 <- ggplot(ssdata1,aes(x=signature,y=Sample,fill=Sig_TCGA)) + xlab("") + ylab("") + geom_tile(color="white",size=0.1) + scale_fill_gradient(low = "gray95", high = "tomato") + geom_text(aes(label=round(Sig_TCGA,2)),angle=45, size=2) +  ggtitle("TCGA") + theme_minimal() +  theme(panel.grid = element_blank(), axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5, size=8))

pp <- ggarrange(p1,p2,ncol=2,nrow=1,common.legend=T,legend="bottom",labels=c("A","B"))
pdf(file="merge.ID.signature.smSNPs.compare.pdf", height=8, width=10)
pp
dev.off()

pdf(file="smsnps.sbs.compare.pdf")
ggplot(data = msdata) + geom_boxplot(mapping = aes(x = signature, y = contribution, color=group),outlier.size = 0.5)  + theme_bw()+ theme(axis.text.y=element_text(size=8), axis.text.x = element_text(size=8), strip.text.y = element_text(size=10), panel.grid = element_blank()) + scale_color_manual(values = colors) +  scale_fill_manual(values = colors) + xlab("signature") + ylab("contribution")  + guides(color = guide_legend( nrow = 2, byrow = TRUE))
dev.off()

