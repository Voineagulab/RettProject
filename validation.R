## load required packages ##################
library(limma)
library(ruv)

## micro-array data #######################
dataExp=read.csv("../Data/dataExpProbes.csv")
infoExp=dataExp[,c(1,2)]
useExp=dataExp[,-c(1,2)]

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")
group_array <- c(rep("control", 6), rep("Rett", 6))

colnames(useExp) <- Samples_array

## Array Nuron Levels Plot ###########################
load("../Data/Astro_array.rda")
barplot(height=1-as.numeric(Astro_array[1,-1]), 
        names.arg=Samples_array,
        ylab = "Proportion of neurons",
        ylim=c(0.48,0.52), col=c(rep("blue", 6), rep("red", 6)),
        xpd=FALSE)
abline(h=0.5, col="grey")

## Array PC plots  ###########################

par(mfrow=c(2,2))
intUQ <- plotMDS(useExp, col=c(rep("BLUE",6),rep("RED",6)), 
                 xlab="PC1", ylab="PC2", 
                 xlim=c(-0.6,1.2), ylim=c(-1.1,0.6),
                 main="UQ",dim.plot = c(1,2), pch=20)
text(intUQ$x, intUQ$y, Samples_array,
     pos=c(1,1,1,1,3,1,3,3,1,2,1,1),
     col=c(rep("BLUE",6),rep("RED",6)))

fit <- eBayes(lmFit(useExp, model.matrix(~group_array)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctl <- !infoExp$TargetID%in%fp_array$TargetID[1:5000]

RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=1)
W2 <- RUVnArray$W
I <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,1]
C <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,-1]
R2 <- useExp - as.matrix(C)%*%as.matrix(t(W2)) - I
intk2 <- plotMDS(R2, col=c(rep("BLUE",6),rep("RED",6)), 
                 xlab="PC1", ylab="PC2", 
                 xlim=c(-0.6,1), ylim=c(-0.6,0.5),
                 main="RUV-2 (k=1)",dim.plot=c(1,2),pch=20)
text(intk2$x, intk2$y, Samples_array,
     pos=c(1,1,1,3,3,1,1,1,1,1,1,1),
     col=c(rep("BLUE",6),rep("RED",6)))

RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=2)
W2 <- RUVnArray$W
I <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,1]
C <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,-1]
R2 <- useExp - as.matrix(C)%*%as.matrix(t(W2)) - I
intk2 <- plotMDS(R2, col=c(rep("BLUE",6),rep("RED",6)), 
                 xlab="PC1", ylab="PC2", 
                 xlim=c(-0.6,1), ylim=c(-0.6,0.5),
                 main="RUV-2 (k=2)",dim.plot=c(1,2),pch=20)
text(intk2$x, intk2$y, Samples_array,
     pos=c(1,1,2,4,3,3,1,1,1,1,1,1),
     col=c(rep("BLUE",6),rep("RED",6)))


RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=3)
W2 <- RUVnArray$W
I <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,1]
C <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,-1]
R2 <- useExp - as.matrix(C)%*%as.matrix(t(W2)) - I
intk2 <- plotMDS(R2, col=c(rep("BLUE",6),rep("RED",6)), 
                 xlab="PC1", ylab="PC2", 
                 xlim=c(-0.6,1), 
                ylim=c(-0.6,0.5),
                 main="RUV-2 (k=3)",dim.plot=c(1,2),pch=20)
text(intk2$x, intk2$y, Samples_array,
     pos=c(1,1,2,4,3,3,1,2,1,2,1,4),
     col=c(rep("BLUE",6),rep("RED",6)))

## RUV2 normalisation with k=2 ##########################
RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=2)
W2 <- RUVnArray$W
I <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,1]
C <- lmFit(useExp, design=model.matrix(~W2))$coefficients[,-1]
R2 <- useExp - as.matrix(C)%*%as.matrix(t(W2)) - I

par(mfrow=c(1,1))
intk2 <- plotMDS(R2, col=c(rep("BLUE",6),rep("RED",6)), 
                 xlab="PC1", ylab="PC2", 
                 xlim=c(-0.6,1), ylim=c(-0.6,0.5),
                 main="RUV-2 (k=2)",dim.plot=c(1,2),pch=20)
text(intk2$x, intk2$y, Samples_array,
     pos=c(1,1,2,4,3,3,1,1,1,1,1,1),
     col=c(rep("BLUE",6),rep("RED",6)))

R2<-cbind(infoExp,R2)
## save(R2, file="../Results/R2.rda")

## Differential expression analysis #################

design_array <- model.matrix(~RUVArray$W+group_array)
fit <- eBayes(lmFit(useExp, design_array))
FDR_array <- topTable(fit, coef=4, number=Inf)
FDR_array <- cbind(infoExp[rownames(FDR_array),], FDR_array[,c(1,2,4,5)])
rownames(FDR_array) <- 1:nrow(FDR_array)

## Adding ProbeID #############################

load("../Data/probe.rda")

FDR_array <- cbind(FDR_array[,1],
                   probe[match(FDR_array$ProbeID, probe$Array_Address_Id),c("Probe_Id")], 
                   FDR_array[,-1])

colnames(FDR_array)[1:3] <- c("GeneName","ProbeID", "AddressID")

FDR_array$GeneName <- as.vector(FDR_array$GeneName)

FDR_array$GeneName[grep("Mar",FDR_array$GeneName)] <- paste("MARCH",substr(FDR_array$GeneName[grep("Mar",FDR_array$GeneName)],1,1),sep="")
FDR_array$GeneName[grep("Sep",FDR_array$GeneName)] <- paste("SEPT",substr(FDR_array$GeneName[grep("Sep",FDR_array$GeneName)],1,1),sep="")

## Adding GeneID ########################

arrayAnnot=read.csv("../Annot/HumanHT-12_V4_0_R2_15002873_B 2.csv")
conv=read.csv("../Annot/conv_refseq_ensembl.csv")
m=match(FDR_array$AddressID, arrayAnnot$Array_Address_Id)
FDR_array$TxID=arrayAnnot$TxID_noV[m]
m=match(FDR_array$TxID, conv$From)
FDR_array$GeneID <- conv$To[m]
FDR_array <- FDR_array[,c("GeneName","GeneID","TxID",
                          "ProbeID","AddressID",
                          "logFC","AveExpr","P.Value","adj.P.Val")]
FDR_array <- FDR_array[!is.na(FDR_array$GeneID),]
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

