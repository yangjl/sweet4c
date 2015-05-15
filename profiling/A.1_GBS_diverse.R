### Jinliang Yang
### May 8th, 2015
### find association for the sweet gene


allzea <- read.csv("~/Documents/Github/SeeDs/data/AllZeaGBSv2.7_publicSamples_metadata20140411.csv")
table(allzea$Project)
table(allzea$GermplasmSet)
ames282 <- subset(allzea, Project %in% "Ames282") 
length(unique(ames282$DNASample)) #283                

ames282$taxa <- paste(ames282$DNASample, ames282$LibraryPrepID, sep=":")
#2010 Ames Lines        AMES Inbreds             Ames282      ApeKI 384-plex 
#2628                 696
write.table(ames282$taxa, "data/Taxa_ames282_288.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
# note this ames lines including 2010 Ames lines and 696 Ames inbred!
# Note: space delimted => %s/\n/ /g in vim