## Load required packages ##########################################
library(edgeR)
library(RUVSeq)
library(EDASeq)

## Load inputs/data #############################################
load("../Results/geneCounts.rda")
load("../Data/Astro_seq.rda")
load("../Annot/geneInfo.rda")

## Two versions of sample naming ###################################

Samples_Seq <- c("C1078", "C1541", "C1571", "R1815", "R4516", "R4852")
Samples_seq <- c("C2", "C3", "C4", "R1", "R2", "R3")

## Normalisation and Filtering ###########################################

geneRPKM <- rpkm(x = geneCounts$counts, gene.length = geneCounts$annotation$Length)
geneCPM <- cpm(geneCounts$counts)
  
threshold = 0.5; nSamples = 3
filtered <- which(rowSums(geneRPKM > threshold) >= nSamples)
counts <- geneCounts$counts[filtered,]
RPKM <- geneRPKM[filtered,]
CPM <- geneCPM[filtered,]
colnames(counts) <- Samples_seq
colnames(RPKM) <- Samples_seq
colnames(CPM) <- Samples_seq

## PC plots ###################################################################
group_seq <- c("control", "control", "control", "Rett", "Rett", "Rett")
input <- DGEList(counts=counts, group=group_seq)
input <- calcNormFactors(input)
plotMDS(input, col=c(rep("blue",3),rep("red",3)),
        xlab="PC1", ylab="PC2", main="TMM Normalisation")
## plotPCA(counts)

## first past differential expression analysis ####################################

input <- estimateCommonDisp(input, verbose=TRUE)
input <- estimateTagwiseDisp(input)
plotBCV(input)
et <- exactTest(input)
firstPass <- topTags(et,n=Inf)
controlGenes <- rownames(firstPass)[-c(1:5000)]

## RUV ##########################################################################

RUV1 <- RUVg(counts, controlGenes, k=1)
N <- RUV1$normalizedCounts
plotMDS(DGEList(RUV1$normalizedCounts), col=c(rep("blue",3),rep("red",3)),
        xlab="PC1", ylab="PC2", main="RUVg (k=1)")

RUV2 <- RUVg(counts, controlGenes, k=2)
plotMDS(DGEList(RUV2$normalizedCounts), col=c(rep("blue",3),rep("red",3)),
        xlab="PC1", ylab="PC2", main="RUVg (k=2)")

## plotPCA(N)
## save(N, file="../Results/N.rda")

## Regions

## Differential expression analysis #############################################
design <- model.matrix(~RUV1$W+group_seq)

input <- DGEList(counts=counts, group=group_seq)
input <- calcNormFactors(input, method="none")

input <- estimateGLMCommonDisp(input, design, verbose=TRUE)
input <- estimateGLMTrendedDisp(input, design)
input <- estimateGLMTagwiseDisp(input, design)
plotBCV(input)

fit <- glmFit(input, design)
et <- glmLRT(fit)

FDR_seq <- data.frame(topTags(et, n=Inf))
FDR_seq <- cbind(rownames(FDR_seq), FDR_seq)
colnames(FDR_seq)[1] <- "genes"
rownames(FDR_seq) <- 1:nrow(FDR_seq)

## Adding Gene Names and Normalised counts ########################################

FDR_seq <- cbind(FDR_seq[,1],
                   geneInfo$gene_name[match(FDR_seq$genes, geneInfo$gene_id)],
                   FDR_seq[,-1])
colnames(FDR_seq)[1:2]<-c("geneID","geneNames")

FDR_seq <- cbind(FDR_seq, N[match(FDR_seq$geneID,rownames(N)),])

## save(FDR_seq, file="../Results/FDR_seq.rda")

nrow(FDR_seq[FDR_seq$FDR<0.05,]) #244 DE genes with FDR<0.05

## MECP2 Expression ##################

## RPKM values
MECP2 <- geneRPKM[rownames(geneRPKM)=="ENSG00000169057"]
barplot(height=as.numeric(MECP2), 
        names.arg=Samples_seq,
        main="MECP2 Expression Level",
        ylab = "RPKM Value",
        ylim=c(7,11), col=c(rep("blue", 3), rep("red", 3)),
        xpd=FALSE)

## RUV Normalised Counts
MECP2 <- FDR_seq[FDR_seq$geneNames=="MECP2",]
##pdf("../Results/MECP2_seq.pdf")
barplot(height=as.numeric(MECP2[-c(1:7)]), 
        names.arg=Samples_seq,
        main="MECP2 Counts",
        ylab = "Normalised Counts",
        ylim=c(800,1500), col=c(rep("blue", 3), rep("red", 3)),
        xpd=FALSE)
##dev.off()

## Regional Analysis #############################

regions <- c("P","T",rep("P",4))

## before RUV
design_regions <- model.matrix(~group_seq+regions)

input <- DGEList(counts=counts)
input <- calcNormFactors(input, method="none")

input <- estimateGLMCommonDisp(input, design_regions, verbose=TRUE)
input <- estimateGLMTrendedDisp(input, design_regions)
input <- estimateGLMTagwiseDisp(input, design_regions)
plotBCV(input)

fit <- glmFit(input, design_regions)
et<-glmLRT(fit,coef=2)

FDR_regions <- data.frame(topTags(et, n=Inf))
FDR_regions <- cbind(rownames(FDR_regions), FDR_regions)
colnames(FDR_regions)[1] <- "genes"
rownames(FDR_regions) <- 1:nrow(FDR_regions)
## 54 regions vs 89 rett/control

## after RUV
design_regions <- model.matrix(~RUV1$W+group_seq+regions)

input <- DGEList(counts=counts)
input <- calcNormFactors(input, method="none")

input <- estimateGLMCommonDisp(input, design_regions, verbose=TRUE)
input <- estimateGLMTrendedDisp(input, design_regions)
input <- estimateGLMTagwiseDisp(input, design_regions)
plotBCV(input)

fit <- glmFit(input, design_regions)

et3 <- glmLRT(fit,coef=3)
FDR_regions3 <- data.frame(topTags(et3, n=Inf))
FDR_regions3 <- cbind(rownames(FDR_regions3), FDR_regions3)
colnames(FDR_regions3)[1] <- "genes"
rownames(FDR_regions3) <- 1:nrow(FDR_regions3)

et4 <- glmLRT(fit,coef=4)
FDR_regions4 <- data.frame(topTags(et4, n=Inf))
FDR_regions4 <- cbind(rownames(FDR_regions4), FDR_regions4)
colnames(FDR_regions4)[1] <- "genes"
rownames(FDR_regions4) <- 1:nrow(FDR_regions4)

## 18 regions vs 223 R-C
## compare with FDR_seq
DEGenes <- FDR_seq[FDR_seq$FDR<0.05,]$geneID
DEGenesp <- FDR_seq[FDR_seq$PValue<0.05,]$geneID
DEGenes2 <- FDR_regions3[FDR_regions3$FDR<0.05,]$genes
DEGenes2p <- FDR_regions3[FDR_regions3$PValue<0.05,]$genes
sig_reg <- FDR_regions4[FDR_regions4$FDR<0.05,]$genes
p_reg <- FDR_regions4[FDR_regions4$PValue<0.05,]$genes

length(DEGenes) #244
length(DEGenes2) #223
length(DEGenesp) #1824
length(DEGenes2p) #1970
length(intersect(DEGenes2,sig_reg)) #5
length(intersect(DEGenes2,p_reg)) #76
length(intersect(DEGenes,DEGenes2)) #135
length(intersect(DEGenesp,DEGenes2p)) #1213
length(intersect(DEGenes,DEGenes2p)) #244
length(intersect(sig_reg,DEGenes)) #2
length(intersect(p_reg,DEGenes)) #18

FDR_seq[FDR_seq$geneID%in%intersect(sig_reg,DEGenes),]$geneNames
FDR_seq[FDR_seq$geneID%in%intersect(p_reg,DEGenes),]$geneNames
FDR_seq[FDR_seq$geneID%in%intersect(DEGenes,DEGenes2),]$geneNames

## Five Sample Analysis ############
group_seq <- c("control", "control", "Rett", "Rett", "Rett")
input <- DGEList(counts=counts[,-2], group=group_seq)
input <- calcNormFactors(input)
plotMDS(input, col=c(rep("blue",2),rep("red",3)),
        xlab="PC1", ylab="PC2", main="TMM Normalisation")
plotPCA(counts[,-2])

input <- estimateCommonDisp(input, verbose=TRUE)
input <- estimateTagwiseDisp(input)
plotBCV(input)
et <- exactTest(input)
firstPass <- topTags(et,n=Inf)
controlGenes <- rownames(firstPass)[-c(1:5000)]

## RUV ##########################################################################

RUV1 <- RUVg(counts[,-2], controlGenes, k=1)
N <- RUV1$normalizedCounts
plotMDS(DGEList(RUV1$normalizedCounts), col=c(rep("blue",2),rep("red",3)),
        xlab="PC1", ylab="PC2", main="RUVg (k=1)")

RUV2 <- RUVg(counts[,-2], controlGenes, k=2)
plotMDS(DGEList(RUV2$normalizedCounts), col=c(rep("blue",2),rep("red",3)),
        xlab="PC1", ylab="PC2", main="RUVg (k=2)")

plotPCA(N)
## save(N, file="../Results/N.rda")

## Differential expression analysis #############################################
design <- model.matrix(~RUV1$W+group_seq)

input <- DGEList(counts=counts[,-2], group=group_seq)
input <- calcNormFactors(input, method="none")

input <- estimateGLMCommonDisp(input, design, verbose=TRUE)
input <- estimateGLMTrendedDisp(input, design)
input <- estimateGLMTagwiseDisp(input, design)
plotBCV(input)

fit <- glmFit(input, design)
et <- glmLRT(fit)

FDR_seq5 <- data.frame(topTags(et, n=Inf))
FDR_seq5 <- cbind(rownames(FDR_seq5), FDR_seq5)
colnames(FDR_seq5)[1] <- "genes"
rownames(FDR_seq5) <- 1:nrow(FDR_seq5)

FDR_seq5 <- cbind(FDR_seq5[,1],
                 geneInfo$gene_name[match(FDR_seq5$genes, geneInfo$gene_id)],
                 FDR_seq5[,-1])
colnames(FDR_seq5)[1:2]<-c("geneID","geneNames")

FDR_seq5 <- cbind(FDR_seq5, N[match(FDR_seq5$geneID,rownames(N)),])

sig5 <- FDR_seq5$geneNames[FDR_seq5$FDR<0.05]
sig6 <- FDR_seq$geneNames[FDR_seq$FDR<0.05]
length(sig5)
length(sig6)
length(intersect(sig5,sig6))


