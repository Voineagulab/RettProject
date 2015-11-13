
## Load packages and inputs ##############################
library(Rsubread)

annotFile="../Annot/genes.gtf"

C1078BAM="../Data/C1078.bam"
C1541BAM="../Data/C1541.bam"
C1571BAM="../Data/C1571.bam"
R1815BAM="../Data/R1815.bam"
R4516BAM="../Data/R4516.bam"
R4852BAM="../Data/R4852.bam"

## Note that C1571BAM has the first entry removed, because featureCounts cannot process the BAM file properly when the first entry is not paired.

## Counting - allow four hours for this step ##########
geneCounts <- featureCounts(c(C1078BAM, C1541BAM, C1571BAM, R1815BAM, R4516BAM, R4852BAM), 
                                  annot.ext=annotFile, isGTFAnnotationFile=TRUE,  
                                  useMetaFeatures=TRUE, GTF.featureType="exon", GTF.attrType="gene_id",
                                  allowMultiOverlap=FALSE, isPairedEnd=TRUE, 
                                  nthreads=8, strandSpecific=2, minMQS=10,
                                  checkFragLength=FALSE, countMultiMappingReads=FALSE, 
                                  requireBothEndsMapped=TRUE,                              
                                  countChimericFragments=FALSE)

## Naming ##################
colnames(geneCounts$counts) <- c("C1078", "C1541", "C1571", "R1815", "R4516", "R4852")
colnames(geneCounts$stat) <- c("Status", "C1078", "C1541", "C1571", "R1815", "R4516", "R4852")

## save(geneCounts, file="../Results/geneCounts.rda") ####
