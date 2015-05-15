### Jinliang Yang
### May 12th, 2015
### 

source("~/Documents/Github/zmSNPtools/Rcodes/mixed_model.R")

    
getpheno <- function(traits="X10KW_Inbred"){
    nathan <- read.csv("data/Nathan_p1_subset.csv")
    a1 <- nrow(nathan)
    nathan$INBRED <- toupper(nathan$INBRED)
    #genotype <- unique(nathan$INBRED)
    #nathan <- read.csv("data/Nathan_p1_subset.csv")
    
    toi <- traits[1]
    nathan <- nathan[nathan[, toi] !=".",]
    nathan[, toi] <- as.numeric(as.character(nathan[, toi]))
    message(sprintf("###>>> total [ %s ], remaining [ %s ] for trait [ %s ]", a1, nrow(nathan), toi))
    out <- mixed_model(data = nathan, model = as.formula(paste(toi, "~ INBRED")), random = ~1 | Env, 
                          trait = toi)
    out$sex <-1
    out <- out[, c(1,3,2)]
    for(i in 2:length(traits)){
        toi <- traits[i]
        nathan <- nathan[nathan[, toi] !=".",]
        nathan[, toi] <- as.numeric(as.character(nathan[, toi]))
        message(sprintf("total [ %s ], remaining [ %s ] for trait [ %s ]", a1, nrow(nathan), toi))
        myblue <- mixed_model(data = nathan, model = as.formula(paste(toi, "~ INBRED")), random = ~1 | Env, 
                              trait = toi)
        out <- merge(out, myblue, by="Genotype", all=TRUE)
    }
    
    nm <- names(out)
    nm <- gsub("_.*", "", nm)
    names(out) <- nm
    return(out)
}

##########
mytraits <- c("X10KW_Inbred","KC_Inbred","CD_Inbred","CL_Inbred","CW_Inbred","TKW_Inbred")
pheno <- getpheno(traits=mytraits)
names(pheno)[1] <- "id"

write.table(pheno, "data/pheno_ames282.txt", row.names=FALSE, quote=FALSE, sep="\t")
