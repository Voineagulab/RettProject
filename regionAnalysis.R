## load required packages ##################
library(limma)
library(ruv)
library(EDASeq)

## micro-array data #######################
load("../Data/useExp.rda")
load("../Data/infoExp.rda")

regionOrder <- c(1,2,5,7,9,11,3,4,6,8,10,12)

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")
group_array <- c(rep("control", 6), rep("Rett", 6))

Samples_Array <- Samples_Array[regionOrder]
Samples_array <- Samples_array[regionOrder]
useExp <- useExp[,regionOrder]
regions <- c(rep("F",6),rep("T",6))

## Array Neuron Levels Plot ###########################
load("../Data/Astro_array.rda")
regionAstro <- Astro_array[1,-1]
regionAstro <- regionAstro[regionOrder]
barplot(height=1-as.numeric(regionAstro), 
        names.arg=Samples_array,
        ylab = "Proportion of neurons",
        ylim=c(0.48,0.52), col=c(rep("blue", 6), rep("red", 6)),
        xpd=FALSE)
abline(h=0.5, col="grey")

plotPCA(as.matrix(useExp),isLog=TRUE)

## RUV2 normalisation with k=2 ##########################
fit <- eBayes(lmFit(useExp, model.matrix(~regions)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctl <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=2)
W <- RUVnArray$W

R <- removeBatchEffect(useExp,covariates=W,design=model.matrix(~regions))
plotPCA(as.matrix(R),isLog=TRUE)
R<-cbind(infoExp,R)

## Frontal vs Temporal #################

design_array <- model.matrix(~W+regions)
fit <- eBayes(lmFit(useExp, design_array))
FDR_array <- topTable(fit, coef=4, number=Inf)
FDR_array <- cbind(infoExp[rownames(FDR_array),], FDR_array[,c(1,2,4,5)])
rownames(FDR_array) <- 1:nrow(FDR_array)

## Frontal vs Temporal - Control ##########################

useExpC <- useExp[,c(1,2,3,7,8,9)]
groupC <- c(rep("F", 3), rep("T", 3))
plotPCA(as.matrix(useExpC),isLog=TRUE)

fit <- eBayes(lmFit(useExpC, model.matrix(~groupC)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctlC <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExpC), X=t(t(c(rep(0,3), rep(1,3)))), ctl=ctlC, k=2)
WC <- RUVnArray$W
RC <- removeBatchEffect(useExpC,covariates=WC,design=model.matrix(~groupC))
plotPCA(as.matrix(RC),isLog=TRUE)
RC<-cbind(infoExp,RC)

design_array <- model.matrix(~WC+groupC)
fit <- eBayes(lmFit(useExpC, design_array))
FDR_C <- topTable(fit, coef=4, number=Inf)
FDR_C <- cbind(infoExp[rownames(FDR_C),], FDR_C[,c(1,2,4,5)])
rownames(FDR_C) <- 1:nrow(FDR_C)
