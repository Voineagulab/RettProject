## load package and data ##########################
library("goseq")
load("../Data/geneInfo.rda")
load("../Results/FDR_seq.rda")
load("../Results/FDR_array.rda")

## Ontology Analysis Functions ############################################

goResults <- function(genes,method="Wallenius") {
  pwf=nullp(genes,"hg19","ensGene")
  GO.wall=goseq(pwf,"hg19","ensGene",method=method)
  go <- cbind(GO.wall, p.adjust(GO.wall$over_represented_pvalue, method="BH"))
  colnames(go)[ncol(go)] <- "FDR"
  go=go[which(go$FDR<0.05),]
  if(nrow(go)>0){
    rownames(go) <- 1:nrow(go)
    genes2go <- getgo(names(genes)[genes==1],"hg19","ensGene")
    go2genes <- goseq:::reversemapping(genes2go)
    go2genes <- go2genes[match(go$category,names(go2genes))]
    
    goGenes <- rep(NA,nrow(go))
    for(i in 1:nrow(go)){
      goGenes[i] <- paste(geneInfo$gene_name[match(go2genes[[i]],geneInfo$gene_id)],collapse=",")
    }
    go<-cbind(go,goGenes)
    colnames(go)[ncol(go)] <- "goGenes"
  }
  return(go)
}


## RNA-Seq ontology analysis #############################

B_seq <- rep(0,nrow(FDR_seq))
names(B_seq) <- FDR_seq$geneID
genes <- B_seq; genes[FDR_seq$FDR<0.05] <- 1
genesDown <- B_seq; genesDown[FDR_seq$FDR<0.05 & FDR_seq$logFC<0] <- 1
genesUp <- B_seq; genesUp[FDR_seq$FDR<0.05 & FDR_seq$logFC>0] <- 1

#pdf("../Results/GO_lengthBias_seq.pdf")
goAll_seq=goResults(genes)
goUp_seq=goResults(genesUp)
goDown_seq=goResults(genesDown)
#dev.off()

## save(goDown_seq,file="../Results/goDown_seq.rda")

## Micro-array ontology analysis ######################################################

B_array_names <- unique(FDR_array$GeneID)
B_array <- rep(0,length(B_array_names))
names(B_array) <- B_array_names

genes <- B_array; genes[names(genes)%in%unique(FDR_array$GeneID[FDR_array$adj.P.Val<0.05])] <- 1
genesDown <- B_array; genesDown[names(genesDown)%in%unique(FDR_array$GeneID[FDR_array$adj.P.Val<0.05 & FDR_array$logFC<0])] <- 1
genesUp <- B_array; genesUp[names(genesDown)%in%unique(FDR_array$GeneID[FDR_array$adj.P.Val<0.05 & FDR_array$logFC>0])] <- 1

## pdf("../Results/GO_lengthBias_array.pdf")
goAll_array=goResults(genes,method="Hypergeometric")
goUp_array=goResults(genesUp,method="Hypergeometric")
goDown_array=goResults(genesDown,method="Hypergeometric")
##dev.off()

## save(goDown_array, file="../Results/goDown_array.rda")

## Comparing the two

goDown_common <- goDown_seq[goDown_seq$category%in%goDown_array$category,c("category","term","ontology")]
goDown_seq_only <- goDown_seq[!goDown_seq$category%in%goDown_array$category,c("category","term","ontology")]
goDown_array_only <- goDown_array[!goDown_array$category%in%goDown_seq$category,c("category","term","ontology")]

