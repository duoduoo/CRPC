setwd("~/Documents/KhuranaLab/EX_0208_prostateRnaCluster")

library(devtools)
library(CMScaller)
library(edgeR)
library(sva)
library(devtools)
library(biomaRt)


################################################################################################################################
# NTP IPM
tmp_org <- read.csv('geneSignature_noPDX_noBatch.csv', header = T)
head(tmp_org)

matrixData <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/IPM_patient_RSEM_combine3MSK.csv', row.names = 1)
emat <- ematAdjust(matrixData, normMethod = 'RLE')
res <- ntp(emat, tmp_org, doPlot=TRUE, nPerm=1000, seed = 42)
View(res)

write.csv(res,'ericaSigGene_ntpPrediction_ipmPatients_noPDX_noBatch.csv', quote = F)

################################################################################################################################
# NTP SU2C
matrixData <- read.csv('su2c_rnaseq_geneNameRevised_noPDX_noBatch_266Sample_new.csv', header = T, row.names = 1)
tmp_org <- read.csv('geneSignature_noPDX_noBatch_hg38_su2c.csv', header = T)
head(tmp_org)

matrixData <- matrixData[rowSums(matrixData) > 0, ]

emat <- ematAdjust(matrixData, center = T, scale = T, normMethod = 'RLE', signalFilt=.1)
res <- ntp(emat, tmp_org, doPlot=TRUE, nPerm=1000, seed = 42)
View(res)
res$prediction <- as.character(res$prediction)
res$prediction[which(res$p.value >= 0.05 | res$FDR >= 0.5)] <- 'unknown'
View(res)
table(res$prediction)
write.csv(res,'ericaSigGene_ntpPrediction_su2cPatients_noPDX_noBatch.csv', quote = F)