### Jinliang Yang
### April 23th, 2015

###
gbs2bed_ames <- function(gbsfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.hmp.txt",
                    outfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames.bed5"){
    
    ### read in GBS file
    #library("data.table")
    ames <- fread(gbsfile, header=TRUE, sep="\t")
    ames <- as.data.frame(ames)
    #message(sprintf("Loaded [ %s ] SNPs and [ %s ] cols for file [%s]!", nrow(gbs), ncol(gbs), gbsfile))
    
    ### change to BED5+ format
    gbs <- ames
    gbs <- gbs[, c(3,4,4,1,2,5, 12:ncol(gbs))]
    names(gbs)[1:6] <- c("chr", "start", "end", "snpid", "alleles", "nchar")
    #nms <- names(gbs)
    #nms <- gsub("\\..*$", "", nms)
    #names(gbs) <- nms
    gbs$start <- gbs$start -1
    message(sprintf("Changed to BED5+ format and start filtering ..."))
    
    ### filter SNPs contain multiple alleles
    gbs$nchar <- nchar(as.character(gbs$alleles))
    subg <- subset(gbs, nchar == 3)
    subg <- subg[, -6]
    #idx <- grep("-", subg$alleles)
    #subg <- subg[-idx,]
    message(sprintf("Remaining [ %s ] sites with two variations!", nrow(subg)))
    
    message(sprintf("Start to IUPAC=>N transforming, recoding and writing ..."))
    
    ###change IUPAC Ambiguity Codes
    #M    A or C    K
    #R    A or G	Y
    #W	A or T	W
    #S	C or G	S
    #Y	C or T	R
    #K	G or T	M
    subg[subg=="M"] <- "N"
    subg[subg=="R"] <- "N"
    subg[subg=="W"] <- "N"
    subg[subg=="S"] <- "N"
    subg[subg=="Y"] <- "N"
    subg[subg=="K"] <- "N"
    
    write.table(subg, outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    message(sprintf("DONE!"))
}





