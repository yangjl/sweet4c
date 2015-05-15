### Jinliang Yang
### May 14th, 2015

####### Plot the missing rate and MAF

frq <- fread("/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.frq", header=TRUE)
frq <- as.data.frame(frq)
frq$chr <- gsub("_.*", "", frq$snpid)
frq$pos <- as.numeric(as.character(gsub(".*_", "", frq$snpid)))

### plot 
pdf("graphs/Figure_ames282_mafmiss.pdf", width=10, height=5)
par(mfrow=c(1,2))
hist(frq$MAF, breaks=50, xlab="Minor Allele Frequency (MAF)", col="bisque2", main="MAF of Diversity Panel (N=282)")
abline(v=0.05, lty=2, lwd=2, col="red")
hist(frq$missing, breaks=50, xlab="Missing Rate (MR)",col="bisque2", main="MR of Diversity Panel (N=282)")
abline(v=0.5, lty=2, lwd=2, col="red")
dev.off()

dim(subset(frq, MAF > 0.05 & missing < 0.5))
### [1] 306206      7
frq1 <- subset(frq, MAF > 0.05 & missing < 0.5)
dim(subset(frq, chr == "S5" & pos > 127464000 & pos < 127468000))
# 3
dim(subset(frq1, chr == "S5" & pos > 127000000 & pos < 128000000))
# 79

