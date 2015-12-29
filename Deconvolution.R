library(EDASeq)
library(edgeR)
library(CellMix)

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
dataExp <- read.csv("../Data/dataExpProbes.csv")
dataExp <- aggregate(dataExp[,-c(1,2)],by=list(dataExp$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])
plotMDS(dataExp)
dec_cahoy=ged(dataExp,ML_cahoy,"ssKL")
dec_array_bRUV <- dec_cahoy@fit@H
#barplot(dec_array_bRUV[1,], ylim=c(0.45, 0.55))
## Part III: Micro-array after RUV ########################################
load("../Results/R1.rda")
dataExp <- aggregate(R1[,-c(1,2)],by=list(R1$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])
plotMDS(dataExp)
dec_cahoy=ged(dataExp,ML_cahoy,"ssKL")
dec_array_aRUV <- dec_cahoy@fit@H
barplot(dec_array_aRUV[1,], ylim=c(0.30, 0.35))
## Boxplots
pdf("BeforeAndAfterRUV_Deconvolution.pdf")
par(mfrow=c(2,1))
boxplot(dec_array_bRUV[2,1:6], dec_array_bRUV[2,7:12], col=c("blue", "red"), ylab="Estimated Proportion of Neurons", main="Before RUV", names=c("Controls", "Rett"))
boxplot(dec_array_aRUV[2,1:6], dec_array_aRUV[2,7:12], col=c("blue", "red"), ylab="Estimated Proportion of Neurons", main="After RUV", names=c("Controls", "Rett"))
dev.off()
