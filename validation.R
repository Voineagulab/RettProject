## load required packages ##################
library(limma)
library(ruv)
library(EDASeq)

## micro-array data #######################
load("../Data/useExp.rda")
load("../Data/infoExp.rda")

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")
group_array <- c(rep("control", 6), rep("Rett", 6))
regions <-c("F","F","T","T","F","T","F","T","F","T","F","T")

## PC Plot ###########################

plotMDS(useExp,col=c(rep("blue",6),rep("red",6)),
        xlab="PC1",ylab="PC2")
## plotPCA(as.matrix(useExp),isLog=TRUE)

## RUV2 normalisation with k=2 ##########################
fit <- eBayes(lmFit(useExp, model.matrix(~group_array)))
fp_array <- topTable(fit, coef=2, number=Inf)
fp_array <- cbind(infoExp[rownames(fp_array),], fp_array[,c(1,2,4,5)])
ctl <- !infoExp$GeneName%in%fp_array$GeneName[1:5000]

RUVnArray <- RUV2(Y=t(useExp), X=t(t(c(rep(0,6), rep(1,6)))), ctl=ctl, k=2)
W <- RUVnArray$W

I <- lmFit(useExp, design=model.matrix(~W))$coefficients[,1]
C <- lmFit(useExp, design=model.matrix(~W))$coefficients[,-1]
R2 <- useExp - as.matrix(C)%*%as.matrix(t(W)) - I
plotMDS(R2, col=c(rep("BLUE",6),rep("RED",6)),
        xlab="PC1", ylab="PC2",xlim=c(-0.6,1), ylim=c(-0.6,0.5),
        main="RUV-2 (k=2)",dim.plot=c(1,2))

R <- removeBatchEffect(useExp,covariates=W,design=model.matrix(~group_array))
## plotMDS(R,col=c(rep("blue",6),rep("red",6)))
## plotPCA(as.matrix(R),isLog=TRUE)
R<-cbind(infoExp,R)
## save(R, file="../Results/R.rda")

## Differential expression analysis #################

design_array <- model.matrix(~W+group_array)
fit <- eBayes(lmFit(useExp, design_array))
FDR_array <- topTable(fit, coef=4, number=Inf)
FDR_array <- cbind(infoExp[rownames(FDR_array),], FDR_array[,c(1,2,4,5)])
rownames(FDR_array) <- 1:nrow(FDR_array)

## Add normalised expression #############
FDR_array <- cbind(FDR_array,R[match(FDR_array$GeneName,R$GeneName),-c(1:5)])
## save(FDR_array, file="../Results/FDR_array.rda") 

## MECP2 Expression ######################

## Raw Expression Level
MECP2<-useExp[infoExp$GeneName=="MECP2",]
barplot(height=as.numeric(MECP2), 
        names.arg=Samples_array,
        main="MECP2 Expression Level",
        ylab = "Expression Level",
        col=c(rep("blue", 6), rep("red", 6)),
        ylim=c(6.5,7.8),
        xpd=FALSE)

## RUV normalised
MECP2 <- FDR_array[FDR_array$GeneName=="MECP2",]
## pdf(file="../Results/MECP2_array.pdf")
barplot(height=as.numeric(MECP2[-c(1:9)]), 
        names.arg=Samples_array,
        main="MECP2 Expression Level",
        ylab = "Normalised Expression Level",
        col=c(rep("blue", 6), rep("red", 6)),
        ylim=c(7,7.8),
        xpd=FALSE)
## dev.off()

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

## Regional Differences ##########

regions <- c("F","F","T","T","F","T","F","T","F","T","F","T")

## before RUV
fit_region <- eBayes(lmFit(useExp, model.matrix(~group_array+regions)))
FDR_region <- topTable(fit_region, coef=3, number=Inf)
FDR_region <- cbind(infoExp[rownames(FDR_region),], FDR_region[,c(1,2,4,5)])

## after RUV
design_region <- model.matrix(~W+group_array+regions)
fit_region <- eBayes(lmFit(useExp, design_region))
FDR_region <- topTable(fit_region, coef=5, number=Inf)
FDR_region <- cbind(infoExp[rownames(FDR_region),], FDR_region[,c(1,2,4,5)])
rownames(FDR_region) <- 1:nrow(FDR_region)

## Regional Differences - Rett Only ##########

regions <- c("F","T","F","T","F","T")

## before RUV
fit_region_bRUV <- eBayes(lmFit(useExp[,7:12], model.matrix(~regions)))
FDR_region_bRUV <- topTable(fit_region_bRUV, coef=2, number=Inf)
FDR_region_bRUV <- cbind(infoExp[rownames(FDR_region_bRUV),], FDR_region_bRUV[,c(1,2,4,5)])

## after RUV
design_region <- model.matrix(~W[7:12,]+regions)
fit_region <- eBayes(lmFit(useExp[,7:12], design_region))
FDR_region <- topTable(fit_region, coef=4, number=Inf)
FDR_region <- cbind(infoExp[rownames(FDR_region),], FDR_region[,c(1,2,4,5)])
rownames(FDR_region) <- 1:nrow(FDR_region)

## Regional Differences - Control Only ##########

regions <- c("F","F","T","T","F","T")

## before RUV
fit_region_bRUV <- eBayes(lmFit(useExp[,1:6], model.matrix(~regions)))
FDR_region_bRUV <- topTable(fit_region_bRUV, coef=2, number=Inf)
FDR_region_bRUV <- cbind(infoExp[rownames(FDR_region_bRUV),], FDR_region_bRUV[,c(1,2,4,5)])

## after RUV
design_region <- model.matrix(~W[1:6,]+regions)
fit_region <- eBayes(lmFit(useExp[,1:6], design_region))
FDR_region <- topTable(fit_region, coef=4, number=Inf)
FDR_region <- cbind(infoExp[rownames(FDR_region),], FDR_region[,c(1,2,4,5)])
rownames(FDR_region) <- 1:nrow(FDR_region)

## Region vs Rett/Control ##############

Region_P <- FDR_region$GeneName[FDR_region$P.Value<0.05]
Region_FDR <- FDR_region$GeneName[FDR_region$adj.P.Val<0.05]
Sig_array <- FDR_array$GeneName[FDR_array$adj.P.Val<0.05]
Sig_seq <- FDR_seq$geneNames[FDR_seq$FDR<0.05]

nrow(FDR_seq) #19470
nrow(FDR_array) #23569
length(Region_P) #1481
length(Region_FDR) #60
length(Sig_array) #318
length(Sig_seq) #244
length(intersect(Region_P,Sig_array)) #38
length(intersect(Region_P,Sig_seq)) #20
length(intersect(Region_FDR,Sig_array)) #0
length(intersect(Region_FDR,Sig_seq))#0

## Differential expression analysis - Frontal ##########################

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

nrow(FDR_F[FDR_F$adj.P.Val<0.05,])
## 28 significant

## Differential expression analysis - Temporal ##########################

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

nrow(FDR_T[FDR_T$adj.P.Val<0.05,])
## 10 significant

## Frontal vs Temporal
sig_F <- FDR_F[FDR_F$adj.P.Val<0.05,]$GeneName
sig_T <- FDR_T[FDR_T$adj.P.Val<0.05,]$GeneName
sig_FT <- intersect(sig_F,sig_T)
## 2 "CSF1R"  "CX3CR1"

