
sdata <- read.table("summary.tcga.cancertype.dbsnp.table.count.txt", header=T,sep="\t")
fisher_test = function(x1,x2,y1,y2){    
    mash = rbind(c(x1,x2),c(y1,y2))
    res = fisher.test(mash)
    return(res$p.value)
}
sdata$fisher_sex_p <- mapply(FUN = fisher_test, sdata$Male, sdata$Female, sdata$T_male, sdata$T_female)
sdata$fisher_age_p <- mapply(FUN = fisher_test, sdata$age1, sdata$age2, sdata$T_age1, sdata$T_age2)
write.table(sdata, file="summary.tcga.cancertype.dbsnp.table.fisher.test.txt", col.names = TRUE, row.names = FALSE, sep="\t",quote=F)


