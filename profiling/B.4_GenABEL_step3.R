### Jinliang Yang
### May 14th, 2015
### plot the simple results and compute the size of effect and heritability

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

load("cache/gwas_res.RData")
res1 <- results(res1)
res1$qval <- p.adjust(res1$P1df, method = "fdr")
res2 <- results(res2)
res2$qval <- p.adjust(res2$P1df, method = "fdr")
res1.eg <- results(res1.eg)
res2.eg <- results(res2.eg)
res1.mm <- results(res1.mm)
res2.mm <- results(res2.mm)

write.table(res1.mm[, 1:11], "data/Table_gwas_10kw.txt", sep="\t", quote=FALSE)
write.table(res2.mm[, 1:11], "data/Table_gwas_tkw.txt", sep="\t", quote=FALSE)

pdf("graphs/Figure_results.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(x=res1.mm$Position, y=-log10(res1.mm$Pc1df), xlab="Chromosome 5", main="Ten Kernel Weight",
     ylab="-log10(q-value)", pch=19, col="cadetblue")
abline(v=127466000, col="blue", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=2)
#S5_128108485 B73=T

plot(x=res2.mm$Position, y=-log10(res2.mm$Pc1df), xlab="Chromosome 5", main="Total Kernel Weight",
     ylab="-log10(q-value)", pch=19, col="cadetblue")
abline(v=127466000, col="blue", lty=2, lwd=2)
abline(h=-log10(0.01), col="red", lty=2, lwd=2)
#S5_127255673 B73=A
dev.off()

######### heritability and effect size #################
geth2 <- function(snpid = "S5_127255673", gm=gm, res1=res1.mm){
    snp <- as.numeric(gtdata(gm[, snpid]))
    #table(snp,exclude=NULL)
    # introduce missing data indicator
    nomissInd <- (!is.na(snp))
    #table(nomissInd)
    gkin <- ibs(gm[nomissInd,], w="freq")
    
    pol1 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    pol2 <- polygenic(X10KW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    #pol1$esth2
    #pol0 <- polygenic(TKW, gm[nomissInd,], kin=gkin, quiet=TRUE)
    #pol0$esth2
    
    ### ===> heritability
    propVarExBySNP1 <- res1[snpid,"chi2.1df"]/res1[snpid,"N"]
    propHerExBySNP1 <- propVarExBySNP1/pol1$esth2
    propVarExBySNP2 <- res1[snpid,"chi2.1df"]/res1[snpid,"N"]
    propHerExBySNP2 <- propVarExBySNP2/pol2$esth2
    res <- data.frame(snp=snpid, trait=c("TKW", "X10KW"), 
                      varp=c(propVarExBySNP1, propVarExBySNP2), 
                      varh2=c(propHerExBySNP1, propHerExBySNP2),
                      h2=c(pol1$esth2, pol2$esth2))
    return(res)
}

h21 <- geth2(snpid = "S5_128108485", gm=gm, res1=res1.mm)
h21 <- h21[-1,]

h22 <- geth2(snpid = "S5_127255673", gm=gm, res1=res2.mm)
h22 <- h22[-2,]
plotbox_1 <- function(){
    snp <- as.numeric(gtdata(gm[, "S5_127255673"]))
    pheno <- phdata(gm)
    df <- merge(pheno, snp, by="row.names")
    df <- subset(df, !is.na(S5_127255673))
    boxplot(TKW ~ as.factor(S5_127255673), data=df, names=c("B73-like","non-B73"))
    print(t.test(subset(df, S5_127255673==0)$TKW, subset(df, S5_127255673==2)$TKW))
}
plotbox_1()

plotbox_2 <- function(){
    snp <- as.numeric(gtdata(gm[, "S5_128108485"]))
    pheno <- phdata(gm)
    df <- merge(pheno, snp, by="row.names")
    df <- subset(df, !is.na(S5_128108485))
    boxplot(X10KW ~ as.factor(S5_128108485), data=df, names=c("B73-like","non-B73"))
}
plotbox_2()