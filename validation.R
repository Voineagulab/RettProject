## load required packages ##################
library(limma)
library(ruv)
library(EDASeq)

## micro-array data #######################
dataExp=read.csv("../Data/dataExpProbes.csv")
infoExp=dataExp[,c(1,2)]
useExp=dataExp[,-c(1,2)]

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")
group_array <- c(rep("control", 6), rep("Rett", 6))

colnames(useExp) <- Samples_array

## Adding ProbeID #############################

load("../Data/probe.rda")
# 
# probe <- data.frame(probe)
# infoExp <- data.frame(infoExp)

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

## Array Neuron Levels Plot ###########################
load("../Data/Astro_array.rda")
barplot(height=1-as.numeric(Astro_array[1,-1]), 
        names.arg=Samples_array,
        ylab = "Proportion of neurons",
        ylim=c(0.48,0.52), col=c(rep("blue", 6), rep("red", 6)),
        xpd=FALSE)
abline(h=0.5, col="grey")

plotPCA(as.matrix(useExp),isLog=TRUE)

## RUV2 normalisation with k=2 ##########################
fit <- eBayes(lmFit(useExp, model.matrix(~group_array)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctl <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=2)
W <- RUVnArray$W

R <- removeBatchEffect(useExp,covariates=W,design=model.matrix(~group_array))
plotPCA(as.matrix(R),isLog=TRUE)
R<-cbind(infoExp,R)
## save(R, file="../Results/R.rda")

## Differential expression analysis #################

design_array <- model.matrix(~W+group_array)
fit <- eBayes(lmFit(useExp, design_array))
FDR_array <- topTable(fit, coef=4, number=Inf)
FDR_array <- cbind(infoExp[rownames(FDR_array),], FDR_array[,c(1,2,4,5)])
rownames(FDR_array) <- 1:nrow(FDR_array)

## save(FDR_array, file="../Results/FDR_array.rda") 

## Calculate validation rate - 74% #################
load("../Results/FDR_seq.rda")

FDR_seq<-FDR_seq[,1:7]
common <- cbind(FDR_seq[,-5],
                FDR_array[match(FDR_seq$geneID, FDR_array$GeneID),-c(1,2)])

colnames(common) <- c("geneID", "geneName",
                      "logFC_seq","logCPM_seq", "pValue_seq", "FDR_seq",
                      "TxID","probeID", "addressID", 
                      "logFC_array", "aveExpr_array","pValue_array", "FDR_array")
common <- common[!is.na(common$probeID),]
rownames(common) <- 1:nrow(common)

DE_seq <- common[common$FDR_seq<0.05,]
validationFDR <- p.adjust(DE_seq$pValue_array,method="BH")
DE_seq <- cbind(DE_seq,validationFDR)
DE_seq_validated <- DE_seq[DE_seq$validationFDR<0.05 & DE_seq$logFC_seq*DE_seq$logFC_array>0, ]
nrow(DE_seq_validated);nrow(DE_seq);nrow(DE_seq_validated)/nrow(DE_seq)

## Differential expression analysis - Frontal ##########################
## 54%
useExpF <- useExp[,c(1,2,5,7,9,11)]
groupF <- c(rep("control", 3), rep("Rett", 3))
fit <- eBayes(lmFit(useExpF, model.matrix(~groupF)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctlF <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExpF), X=t(t(c(rep(0,3), rep(1,3)))), ctl=ctlF, k=2)
WF <- RUVnArray$W
RF <- removeBatchEffect(useExpF,covariates=WF,design=model.matrix(~groupF))
plotPCA(as.matrix(RF),isLog=TRUE)
RF<-cbind(infoExp,RF)

design_array <- model.matrix(~WF+groupF)
fit <- eBayes(lmFit(useExpF, design_array))
FDR_F <- topTable(fit, coef=4, number=Inf)
FDR_F <- cbind(infoExp[rownames(FDR_F),], FDR_F[,c(1,2,4,5)])
rownames(FDR_F) <- 1:nrow(FDR_F)

common <- cbind(FDR_seq[,-5],
                FDR_F[match(FDR_seq$geneID, FDR_F$GeneID),-c(1,2)])

colnames(common) <- c("geneID", "geneName",
                      "logFC_seq","logCPM_seq", "pValue_seq", "FDR_seq",
                      "TxID","probeID", "addressID", 
                      "logFC_array", "aveExpr_array","pValue_array", "FDR_array")
common <- common[!is.na(common$probeID),]
rownames(common) <- 1:nrow(common)

DE_seq <- common[common$FDR_seq<0.05,]
validationFDR <- p.adjust(DE_seq$pValue_array,method="BH")
DE_seq <- cbind(DE_seq,validationFDR)
DE_seq_VF <- DE_seq[DE_seq$validationFDR<0.05 & DE_seq$logFC_seq*DE_seq$logFC_array>0, ]
nrow(DE_seq_VF);nrow(DE_seq);nrow(DE_seq_VF)/nrow(DE_seq)

## Differential expression analysis - Temporal ##########################
## 67%
useExpT <- useExp[,c(3,4,6,8,10,12)]
fit <- eBayes(lmFit(useExpT, model.matrix(~groupF)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctlT <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExpT), X=t(t(c(rep(0,3), rep(1,3)))), ctl=ctlT, k=2)
WT <- RUVnArray$W
RT <- removeBatchEffect(useExpT,covariates=WT,design=model.matrix(~groupF))
plotPCA(as.matrix(RT),isLog=TRUE)
RT<-cbind(infoExp,RT)

design_array <- model.matrix(~WT+groupF)
fit <- eBayes(lmFit(useExpT, design_array))
FDR_T <- topTable(fit, coef=4, number=Inf)
FDR_T <- cbind(infoExp[rownames(FDR_T),], FDR_T[,c(1,2,4,5)])
rownames(FDR_T) <- 1:nrow(FDR_T)

common <- cbind(FDR_seq[,-5],
                FDR_T[match(FDR_seq$geneID, FDR_T$GeneID),-c(1,2)])

colnames(common) <- c("geneID", "geneName",
                      "logFC_seq","logCPM_seq", "pValue_seq", "FDR_seq",
                      "TxID","probeID", "addressID", 
                      "logFC_array", "aveExpr_array","pValue_array", "FDR_array")
common <- common[!is.na(common$probeID),]
rownames(common) <- 1:nrow(common)

DE_seq <- common[common$FDR_seq<0.05,]
validationFDR <- p.adjust(DE_seq$pValue_array,method="BH")
DE_seq <- cbind(DE_seq,validationFDR)
DE_seq_VT <- DE_seq[DE_seq$validationFDR<0.05 & DE_seq$logFC_seq*DE_seq$logFC_array>0, ]
nrow(DE_seq_VT);nrow(DE_seq);nrow(DE_seq_VT)/nrow(DE_seq)

## Frontal vs Temporal ############
Z<-intersect(DE_seq_VF$geneName,DE_seq_VT$geneName)
q<-length(Z)
m<-nrow(DE_seq_VF)
n<-nrow(DE_seq)-m
k<-nrow(DE_seq_VT)
phyper(q,m,n,k,lower.tail = FALSE)

## complement cascade #######
CC <- c("C1QA","C1QB","C1QC",
        "C3","TGFBR2","CX3CR1","TYROBP")
CC%in%Z
CC%in%DE_seq$geneName
CC%in%DE_seq_VT$geneName
CC%in%DE_seq_VF$geneName

## More comparisons ############
A <- unique(union(DE_seq_VT$geneName,DE_seq_VF$geneName))
length(A)/nrow(DE_seq)
## 72%
B<-unique(union(A,DE_seq_validated))
length(B)/nrow(DE_seq)
## 80%
C<-setdiff(DE_seq_validated$geneName,union(DE_seq_VT$geneName,DE_seq_VF$geneName))
length(C)
D<-setdiff(union(DE_seq_VT$geneName,DE_seq_VF$geneName),DE_seq_validated$geneName)
length(D)
setdiff(intersect(DE_seq_VT$geneName,DE_seq_VF$geneName),DE_seq_validated$geneName)
