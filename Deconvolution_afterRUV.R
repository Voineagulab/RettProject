## Part I: Micro-array Data ########################################

library(CellMix)

## RUV normalised micro-array data ################################################
load("../Results/R1.rda")

dataExp <- aggregate(R1[,-c(1,2)],by=list(R1$TargetID),FUN=mean)
rownames(dataExp)<-dataExp$Group.1
dataExp<-as.matrix(dataExp[,-1])

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
sig <- data.frame(sig)

## RNA-Seq Data #############################################
load("../Results/FDR_seq.rda")
counts<-FDR_seq[,8:13]
geneTPM <- cpm(counts)

TPM<-data.frame(geneTPM)
TPM <- cbind(FDR_seq$geneNames,TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]

common <- intersect(rownames(TPM),rownames(sig))
TPM<-TPM[rownames(TPM)%in%common,]

sig <- data.frame(sig[match(rownames(TPM),rownames(sig)),])

sig <- sig[sig$Neurons>100 & sig$Astrocytes>100,]
sig <- sig[sig$Neurons<1000&
             sig$Astrocytes<1000,]
TPM <- TPM[match(rownames(sig),rownames(TPM)),]
## Decovlution ###################################
output <- DeconRNASeq(datasets=TPM, signatures=sig)
np_afterRUV <- output$out.all[,1]

## Cahoy #######################################
load("../Results/FDR_seq.rda")
counts<-FDR_seq[,8:13]
geneTPM <- cpm(counts)

TPM<-data.frame(geneTPM)
TPM <- cbind(FDR_seq$geneNames,TPM)
colnames(TPM)[1]<-"geneName"
TPM <- aggregate(TPM[,-1],by=list(TPM$geneName),FUN=mean)
rownames(TPM) <- TPM$Group.1
TPM <-TPM[,-1]
TPM<-as.matrix(TPM)
  
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
