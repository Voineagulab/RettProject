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
## Part IV: RNA-Seq before RUV #######################################

load("../Results/geneCounts.rda")
load("../Annot/geneInfo.rda")

colnames(geneCounts$counts) <- c("C2","C3","C4","R1","R2","R3")
geneRPKM <- rpkm(x = geneCounts$counts, gene.length = geneCounts$annotation$Length)
geneTPM <- cpm(geneCounts$counts)

threshold = 0.5; nSamples = 3
filtered <- which(rowSums(geneRPKM > threshold) >= nSamples)
geneCounts <- geneCounts$counts[filtered,]
plotPCA(geneCounts)

TPM <- geneTPM[filtered,]
TPM<-data.frame(TPM)
TPM <- cbind(geneInfo$gene_name[match(rownames(TPM),geneInfo$gene_id)],TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]
TPM <- as.matrix(TPM)

dec_cahoy=ged(TPM,ML_cahoy,"ssKL")
dec_seq_bRUV <- dec_cahoy@fit@H

## Part V: RNA-Seq after RUV #######################################
load("../Results/FDR_seq.rda")
geneCounts<-FDR_seq[,8:13]
plotPCA(as.matrix(geneCounts))
geneTPM <- cpm(geneCounts)

TPM<-data.frame(geneTPM)
TPM <- cbind(FDR_seq$geneNames,TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]
TPM<-as.matrix(TPM)

dec_cahoy=ged(TPM,ML_cahoy,"ssKL")
dec_seq_aRUV <- dec_cahoy@fit@H
