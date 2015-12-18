## Part I: Micro-array Data ########################################

library(CellMix)

## micro-array data ################################################
dataExp <- read.csv("../Data/dataExpProbes.csv")
dataExp <- aggregate(dataExp[,-c(1,2)],by=list(dataExp$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])

## Torkanmani ####################################

tork=read.csv("../Data/Torkamani.csv", header = T)
asts=toupper(tork$Astrocytes[!tork$Astrocytes==""])
neus=toupper(tork$Neurons[!tork$Neurons==""])
common <- intersect(asts,neus)
asts <- unique(setdiff(asts,common))
neus <- unique(setdiff(neus,common))

txtdesc <- function(x) textConnection(paste(x, collapse="\n"))
ML_tork=MarkerList(file=txtdesc(c(paste(asts,"Astrocytes"),paste(neus,"Neurons"))))

dec_tork=ged(dataExp,ML_tork,"ssKL")
dec_tork_ssKL <- dec_tork@fit@H

## Cahoy #########################################################
cahoy_ast=read.csv("../Data/Cahoy_S4_ast.csv", header = T)
cahoy_neu=read.csv("../Data/Cahoy_S6_neu.csv", header = T)

asts <- toupper(cahoy_ast$Gene.Name)
neus <- toupper(cahoy_neu$Gene.Name)
common <- intersect(asts,neus)
asts <- unique(setdiff(asts,common))
neus <- unique(setdiff(neus,common))

txtdesc <- function(x) textConnection(paste(x, collapse="\n"))
ML_cahoy=MarkerList(file=txtdesc(c(paste(asts,"Astrocytes"),paste(neus,"Neurons"))))

dec_cahoy=ged(2^dataExp,ML_cahoy,"ssKL")
dec_cahoy_ssKL <- dec_cahoy@fit@H

## Part II: RNA-Seq Data #######################################

library(DeconRNASeq)
library(edgeR)

## signatures ##################
signatures <- read.csv("../Data/hg19.cage_peak_tpm_ann.osc_pcells_brainExpressed_0.5TPMx2_aggregatedGene.csv")
SampleInfo <- read.csv("../Data/Fantom5SamplesAnnotated.csv")
n<- grep("Neurons",SampleInfo$Sample)
a<- grep("Astr",SampleInfo$Sample)

nL <- SampleInfo$Library[n[1:3]]
nS <- signatures[,colnames(signatures)%in%nL]
aL <- SampleInfo$Library[a[4:6]]
aS <- signatures[,colnames(signatures)%in%aL]

sig <- cbind(rowMeans(nS),rowMeans(aS))
rownames(sig) <- signatures$X
colnames(sig) <- c("Neurons","Astrocytes")

## RNA-Seq Data ##############

load("../Results/geneCounts.rda")
load("../Annot/geneInfo.rda")

geneRPKM <- rpkm(x = geneCounts$counts, gene.length = geneCounts$annotation$Length)
geneTPM <- cpm(geneCounts$counts)

threshold = 0.5; nSamples = 3
filtered <- which(rowSums(geneRPKM > threshold) >= nSamples)
TPM <- geneTPM[filtered,]

TPM<-data.frame(TPM)
TPM <- cbind(geneInfo$gene_name[match(rownames(TPM),geneInfo$gene_id)],TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]

common <- intersect(rownames(TPM),rownames(sig))
TPM<-TPM[rownames(TPM)%in%common,]

sig <- data.frame(sig[match(rownames(TPM),rownames(sig)),])

sig <- sig[sig$Neurons>50 & sig$Astrocytes>50,]
sig <- sig[sig$Neurons<1000&
             sig$Astrocytes<1000,]
TPM <- TPM[match(rownames(sig),rownames(TPM)),]

## Decovlution ###################################
output <- DeconRNASeq(datasets=TPM, signatures=sig)
np <- output$out.all[,1]

## Cahoy #################################

load("../Results/geneCounts.rda")
load("../Annot/geneInfo.rda")

geneRPKM <- rpkm(x = geneCounts$counts, gene.length = geneCounts$annotation$Length)
geneTPM <- cpm(geneCounts$counts)

threshold = 0.5; nSamples = 3
filtered <- which(rowSums(geneRPKM > threshold) >= nSamples)
TPM <- geneTPM[filtered,]

TPM<-data.frame(TPM)
TPM <- cbind(geneInfo$gene_name[match(rownames(TPM),geneInfo$gene_id)],TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]
TPM <- as.matrix(TPM)

cahoy_ast=read.csv("../Data/Cahoy_S4_ast.csv", header = T)
cahoy_neu=read.csv("../Data/Cahoy_S6_neu.csv", header = T)

asts <- toupper(cahoy_ast$Gene.Name)
neus <- toupper(cahoy_neu$Gene.Name)
common <- intersect(asts,neus)
asts <- unique(setdiff(asts,common))
neus <- unique(setdiff(neus,common))

txtdesc <- function(x) textConnection(paste(x, collapse="\n"))
ML_cahoy=MarkerList(file=txtdesc(c(paste(asts,"Astrocytes"),paste(neus,"Neurons"))))

dec_cahoy=ged(TPM,ML_cahoy,"ssKL")
dec_cahoy_ssKL <- dec_cahoy@fit@H


## Cahoy #################################

load("../Results/geneCounts.rda")
load("../Annot/geneInfo.rda")

geneRPKM <- rpkm(x = geneCounts$counts, gene.length = geneCounts$annotation$Length)
geneTPM <- cpm(geneCounts$counts)

threshold = 0.5; nSamples = 3
filtered <- which(rowSums(geneRPKM > threshold) >= nSamples)
counts <- geneCounts$counts[filtered,]
TPM <- 2^(voom(counts)$E)

TPM<-data.frame(TPM)
TPM <- cbind(geneInfo$gene_name[match(rownames(TPM),geneInfo$gene_id)],TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]
TPM <- as.matrix(TPM)

cahoy_ast=read.csv("../Data/Cahoy_S4_ast.csv", header = T)
cahoy_neu=read.csv("../Data/Cahoy_S6_neu.csv", header = T)

asts <- toupper(cahoy_ast$Gene.Name)
neus <- toupper(cahoy_neu$Gene.Name)
common <- intersect(asts,neus)
asts <- unique(setdiff(asts,common))
neus <- unique(setdiff(neus,common))

txtdesc <- function(x) textConnection(paste(x, collapse="\n"))
ML_cahoy=MarkerList(file=txtdesc(c(paste(asts,"Astrocytes"),paste(neus,"Neurons"))))

dec_cahoy=ged(TPM,ML_cahoy,"ssKL")
dec_cahoy_ssKL <- dec_cahoy@fit@H



