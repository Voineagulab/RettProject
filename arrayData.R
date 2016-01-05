
## micro-array data #######################
dataExp=read.csv("../Data/dataExpProbes.csv")
infoExp=dataExp[,c(1,2)]
useExp=dataExp[,-c(1,2)]

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")
group_array <- c(rep("control", 6), rep("Rett", 6))
regions <-c("F","F","T","T","F","T","F","T","F","T","F","T")

colnames(useExp) <- Samples_array

## Adding ProbeID #############################

load("../Data/probe.rda")

infoExp <- data.frame(infoExp[,1],
                 probe[match(infoExp$ProbeID, probe$Array_Address_Id),c("Probe_Id")], 
                 infoExp[,-1])

colnames(infoExp)[1:3] <- c("GeneName","ProbeID", "AddressID")

## Correcting gene names ############

infoExp$GeneName <- as.vector(infoExp$GeneName)

infoExp$GeneName[grep("Mar",infoExp$GeneName)] <- paste("MARCH",substr(infoExp$GeneName[grep("Mar",infoExp$GeneName)],1,1),sep="")
infoExp$GeneName[grep("Sep",infoExp$GeneName)] <- paste("SEPT",substr(infoExp$GeneName[grep("Sep",infoExp$GeneName)],1,1),sep="")

## Adding GeneID ########################

arrayAnnot=read.csv("../Annot/HumanHT-12_V4_0_R2_15002873_B 2.csv")
conv=read.csv("../Annot/conv_refseq_ensembl.csv")
m=match(infoExp$AddressID, arrayAnnot$Array_Address_Id)
infoExp$TxID=arrayAnnot$TxID_noV[m]
m=match(infoExp$TxID, conv$From)
infoExp$GeneID <- conv$To[m]
infoExp <- infoExp[,c("GeneName","GeneID","TxID","ProbeID","AddressID")]

## save(infoExp, file="../Data/infoExp.rda")
## save(useExp, file="../Data/useExp.rda")