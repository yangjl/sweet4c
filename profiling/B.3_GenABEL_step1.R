### Jinliang Yang
### May 12th, 2015
### find the common set of the genotype and phenotype data

#### Germplasm GWAS using GenABEL
commonGenoPheno <- function(){
    #### ==>
    pheno <- read.table("data/pheno_ames282.txt", header=TRUE)
    
    #### ==>
    geno <- fread("largedata/4.sweet/ZeaGBSv27_Ames282.illumina", sep="\t", header=TRUE)
    geno <- as.data.frame(geno)
    
    names(geno) <- toupper(names(geno))
    idx <- names(geno)[names(geno) %in% pheno$id]
    message(sprintf("###>>> common ids [ %s ]", length(idx)))
    geno <- geno[, c("NAME", "CHR", "POS", idx)]
    names(geno)[1:3] <- c("Name", "Chr", "Pos")
    write.table(geno, "largedata/4.sweet/ZeaGBSv27_Ames264.illumina", sep="\t", row.names=FALSE, quote=FALSE)
    
}

library("data.table", lib="~/bin/Rlib/")
library("GenABEL.data", lib="~/bin/Rlib/")
library("GenABEL", lib="~/bin/Rlib/")

commonGenoPheno()
###>>> common ids [ 264 ]
convert.snp.illumina(infile="largedata/4.sweet/ZeaGBSv27_Ames264.illumina", 
                     outfile="largedata/4.sweet/geno_ames264.raw", strand = "+", bcast = 10000000)



