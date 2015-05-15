### Jinliang Yang
### purpose: transform ZeaGBSv2.7 Ames inbred => bed+ format

######=======>
library(data.table, lib = "~/bin/Rlib")

source("lib/gbs2bed_ames282.R")
gbs2bed_ames(gbsfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.hmp.txt",
             outfile="/group/jrigrp4/AllZeaGBSv2.7impV5/ZeaGBSv27_Ames282.bed5")

### log
#>>> Read 955690 rows and 299 (of 299) columns from 0.555 GB file in 00:00:34
#>>> Changed to BED5+ format and start filtering ...
#>>> Remaining [ 509572 ] sites with two variations!
#>>> Start to IUPAC=>N transforming, recoding and writing ...

