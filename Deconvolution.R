## Set-up ##########################
library(EDASeq)
library(edgeR)
library(CellMix)

Samples_Array <-  c("C55F","C1078F","C1078T","C1541T","C1571F","C1571T","R1815F","R1815T","R4516F","R4516T","R4852F","R4852T")
Samples_array <- c("C1-F", "C2-F", "C2-T", "C3-T", "C4-F", "C4-T", "R1-F", "R1-T", "R2-F", "R2-T", "R3-F", "R3-T")

## Part I: Markers - Cahoy ##############################################
cahoy_ast=read.csv("../Data/Cahoy_S4_ast.csv", header = T)
cahoy_neu=read.csv("../Data/Cahoy_S6_neu.csv", header = T)

asts <- toupper(cahoy_ast$Gene.Name)
neus <- toupper(cahoy_neu$Gene.Name)
common <- intersect(asts,neus)
asts <- unique(setdiff(asts,common))
neus <- unique(setdiff(neus,common))

txtdesc <- function(x) textConnection(paste(x, collapse="\n"))
ML_cahoy=MarkerList(file=txtdesc(c(paste(asts,"Astrocytes"),paste(neus,"Neurons"))))

## Part II: Micro-array before RUV ########################################
pdf("../Results/Plots_DeconvArray.pdf", height=10, width = 10)
## log scale
dataExp <- read.csv("../Data/dataExpProbes.csv")
colnames(dataExp)[-c(1:2)] <- Samples_array
dataExp <- aggregate(dataExp[,-c(1,2)],by=list(dataExp$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])
dec_cahoy=ged(dataExp,ML_cahoy,"ssKL")
dec_array_bRUV <- dec_cahoy@fit@H
barplot(height=as.numeric(dec_array_bRUV[2,]), 
        names.arg=colnames(dec_array_bRUV),
        ylab = "Neuron Percentage",
        ylim=c(0.46,0.51), col=c(rep("blue", 6), rep("red", 6)),
        xpd=FALSE,  main="Deconv before RUV,log2 transformed data")

## Non-log scale
dataexp <- read.csv("/Users/rna/Google Drive/DATA/AMELIA ASSAREH/DECONVOLUTION/Data/dataExpProbes_noLog.csv", head = T)
colnames(dataexp)[-c(1:2)] <- Samples_array
dataexp <- aggregate(dataexp[,-c(1,2)],by=list(dataexp$TargetID),FUN=mean)
rownames(dataexp)<-dataexp$Group.1
dataexp<-as.matrix(dataexp[,-1])
dec_cahoy=ged(dataexp,ML_cahoy,"ssKL")
dec_array_bRUV <- dec_cahoy@fit@H
barplot(height=as.numeric(dec_array_bRUV[2,]), 
        names.arg=colnames(dec_array_bRUV),
        ylab = "Neuron Percentage",
        ylim=c(0.40,0.60), col=c(rep("blue", 6), rep("red", 6)),
        xpd=FALSE,  main="Deconv before RUV, NON-log2 transformed data")

## Part III: Micro-array after RUV ########################################
load("../Results/R1.rda")

dataExp <- aggregate(R1[,-c(1,2)],by=list(R1$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])
dec_cahoy=ged(2^dataExp,ML_cahoy,"ssKL")
dec_array_aRUV <- dec_cahoy@fit@H
barplot(dec_array_aRUV[2,], col=c(rep("blue",6),rep("red",6)),ylim=c(0.59, 0.66),xpd=FALSE, main="Deconv after RUV,NON-log2 transformed data")
## Boxplots
par(mfrow=c(2,1))
boxplot(dec_array_bRUV[2,1:6], dec_array_bRUV[2,7:12], col=c("blue", "red"), ylab="Estimated Proportion of Neurons", main="Before RUV", names=c("Controls", "Rett"))
boxplot(dec_array_aRUV[2,1:6], dec_array_aRUV[2,7:12], col=c("blue", "red"), ylab="Estimated Proportion of Neurons", main="After RUV", names=c("Controls", "Rett"))
dev.off()
