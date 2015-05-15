### Jinliang Yang
### May 12th, 2015
### conduct MLM GWAS

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

### load the data
gm <- load.gwaa.data(phe="data/pheno_ames282.txt", gen="largedata/4.sweet/geno_ames264.raw", force=T)
head(gm@phdata)
#gm@gtdata

################# QC1 #########################
qc1 <- check.marker(gm, p.level=0, callrate=0.5, ibs.mrk=0, maf=0.05, perid.call=0.5)
# 0 people excluded because too high autosomal heterozygosity (FDR <1%)
# In total, 304009 (100%) markers passed all criteria
# In total, 264 (100%) people passed all criteria

### using the simple model GWAS
frq <- fread("largedata/4.sweet/ZeaGBSv27_Ames264.illumina", header=TRUE)
frq <- as.data.frame(frq)
frq <- frq[, 1:3]
### gene located on 127466k region.
chr5snp <- subset(frq, Chr == 5 &Pos > 122000000 & Pos < 132000000)$Name
#349
gm.qc1 <- gm[qc1$idok, chr5snp]

############### simple model ##########################
par(mfrow=c(1,2))
res1 <- qtscore(X10KW, data=gm.qc1, trait = "gaussian" )
res2 <- qtscore(TKW, data=gm.qc1, trait = "gaussian")

plot(res1)
abline(v=127466000)
plot(res2)
abline(v=127466000)

############### kinship ##########################
gkin <- ibs(gm, weight="freq")
res1.eg <- egscore(X10KW, data=gm.qc1, kinship.matrix=gkin)
res2.eg <- egscore(TKW, data=gm.qc1, kinship.matrix=gkin)

plot(res1.eg)
abline(v=127466000)
plot(res2.eg)
abline(v=127466000)

############### MLM #################
h2a1 <- polygenic(X10KW, data=gm, kin=gkin)
res1.mm <- mmscore(h2a1, data=gm, snpsubset=chr5snp)

h2a2 <- polygenic(TKW, data=gm, kin=gkin)
res2.mm <- mmscore(h2a2, data=gm, snpsubset=chr5snp)

plot(res1.mm, main="10 kernel weight", pch=16, col="cadetblue")
abline(v=127466000)
plot(res2.mm, main="total kernel weight", pch=16, col="cadetblue")
abline(v=127466000)
save(file="cache/gwas_res.RData", list=c("gm", "res1", "res2", "res1.eg", "res2.eg", "res1.mm", "res2.mm"))