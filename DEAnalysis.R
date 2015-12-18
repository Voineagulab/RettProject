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

## Neuron Levels Plot ####################################################### 

barplot(height=1-as.numeric(Astro_seq), 
        names.arg=Samples_seq,
        ylab = "Proportion of Neurons",
        ylim=c(0.48,0.52), col=c(rep("red", 3), rep("blue", 3)),
        xpd=FALSE)
abline(h=0.5, col="grey")

## PC plots ###################################################################
group_seq <- c("control", "control", "control", "Rett", "Rett", "Rett")
input <- DGEList(counts=counts, group=group_seq)
input <- calcNormFactors(input)
plotMDS(input, col=c(rep("blue",3),rep("red",3)),
        xlab="PC1", ylab="PC2", main="TMM Normalisation")
plotPCA(counts)

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

plotPCA(N)
## save(N, file="../Results/N.rda")
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

