### Jinliang
### May 12th, 2015

#source("~/Documents/Github/zmSNPtools/Rcodes/dsnp2GenABEL.R")
library("data.table", lib="~/bin/Rlib/")
#library("GenABEL.data", lib="~/bin/Rlib/")
#library("GenABEL", lib="~/bin/Rlib/")

bed2illumina <- function(mafcutoff=0.05, missingrate=0.5){
    #######==> GBS data for GenABEL
    gbs <- fread("/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.bed5")
    gbs <- as.data.frame(gbs)
    
    frq <- fread("/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.frq", header=TRUE)
    frq <- as.data.frame(frq)
    frq$chr <- gsub("_.*", "", frq$snpid)
    #frq$pos <- as.numeric(as.character(gsub(".*_", "", frq$snpid)))
    
    flt <- subset(frq, chr != "S0" & MAF > mafcutoff & missing < missingrate)
    #dim(subset(flt, chr == "S5" & pos > 127464000 & pos < 127468000))
    # 3
    #dim(subset(flt, chr == "S5" & pos > 127000000 & pos < 128000000))
    # 135
    gbs1 <- subset(gbs, snpid %in% flt$snpid)
    message(sprintf("Input [ %s ] GBS data, after filtering, [ %s ] remaining", nrow(gbs), nrow(flt)))
    nms <- names(gbs1)
    nms2 <- gsub(":.*", "", nms)
    names(gbs1) <- nms2
    
    gbs1 <- gbs1[, c(-2, -5)]
    gbs1 <- gbs1[, c(3,1, 2, 4:ncol(gbs1))]
    names(gbs1)[1:3] <- c("Name", "Chr", "Pos")
    
    gbs1[gbs1 == "A"] <- "AA"
    gbs1[gbs1 == "T"] <- "TT"
    gbs1[gbs1 == "C"] <- "CC"
    gbs1[gbs1 == "G"] <- "GG"
    gbs1[gbs1 == "-"] <- "--"
    gbs1[gbs1 == "N"] <- "00"
    
    message("start to writing ...")
    write.table(gbs1, "largedata/4.sweet/ZeaGBSv27_Ames282.illumina", sep="\t", row.names=FALSE, quote=FALSE)
}

bed2illumina(mafcutoff=0.05, missingrate=0.5)
###>>> Read 509572 rows and 293 (of 293) columns from 0.291 GB file in 00:00:14
###>>> Input [ 509572 ] GBS data, after filtering, [ 306190 ] remaining

