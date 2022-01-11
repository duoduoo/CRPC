setwd("~/Documents/KhuranaLab/EX_0209_crpcFigures")

library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#### Figure 1D #############################################################
## feature distribution of peaks

## feature distribution idr peaks w PDX but droping the last 3 samples
peaks <- dir("../EX_0201_prostateSubtypePredict/IDR_peaks/", pattern = "*.narrowPeak$", 
             full.names = TRUE)
organoid_peaks_all <- as.list(peaks)
organoid_peaks_all
test <- unlist(organoid_peaks_all)
test2 <- as.vector(sapply(test, function(x) unlist(strsplit(x, "IDR_peaks//"))[2]))
test3 <- as.vector(sapply(test2, function(x) unlist(strsplit(x, ".idr"))[1]))
test4 <- as.vector(sapply(test3, function(x) unlist(strsplit(x, "_"))[1]))
test4
test4[2] <- 'MSKPCa18'
test4[3] <- 'MSKPCa19'
test4[c(40,41,42)] <- c('MSKPCa20','MSKPCa22','MSKPCa24')
test4[4] <- 'C4-2'
test4[36:39] <- c("WCM1078","WCM1262","WCM154", "WCM155")
test4
names(organoid_peaks_all) <- test4

sampleReordered <- organoid_peaks_all[c("MSKPCa1", "MSKPCa2", "MSKPCa3", "MSKPCa8", "MSKPCa9", "MSKPCa10", "MSKPCa11", "MSKPCa12", "MSKPCa13", "MSKPCa14", "MSKPCa15", "MSKPCa16", "MSKPCa17", "MSKPCa18", "MSKPCa19", "MSKPCa20", "MSKPCa22", "MSKPCa24", "WCM154", "WCM155", "WCM1078", "WCM1262", "LNCaP", "VCaP", "22Rv1", "C4-2", "PC3", "DU145", "H660","PDX01a","PDX09a","PDX10a","PDX16a","PDX16b","PDX34a", "MSKEF1", "PARCB1", "PARCB3", "PARCB6", "PARCB8")]

peakAnnoList <- lapply(sampleReordered, annotatePeak, 
                       tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Hs.eg.db", TxDb=txdb)
saveRDS(peakAnnoList, 'figure1d_peakAnnoList_wPDX_drop3pdx.RDS')


pdf("figure1/figure1d_wPDX_drop3pdx.pdf", width=6, height=6, useDingbats=F)
plotAnnoBar(peakAnnoList,
            title = "Feature Distribution")
dev.off()

peakAnnoList <- readRDS('figure1d_peakAnnoList_wPDX_drop3pdx.RDS')
for(x in 1:40){
  print(peakAnnoList[x])
}

peakAnnoList[1][1]
test <- plotAnnoBar(peakAnnoList,
                    title = "Feature Distribution")
write.csv(test$data, 'tables/peakAnnotate.csv')

#################################################################
## 

#################################################################
## feature distribution idr peaks wo PDX
peaks <- dir("../EX_0201_prostateSubtypePredict/IDR_peaks/", pattern = "*.narrowPeak$", 
             full.names = TRUE)
organoid_peaks_all <- as.list(peaks)
organoid_peaks_all
test <- unlist(organoid_peaks_all)
test2 <- as.vector(sapply(test, function(x) unlist(strsplit(x, "IDR_peaks//"))[2]))
test3 <- as.vector(sapply(test2, function(x) unlist(strsplit(x, ".idr"))[1]))
test4 <- as.vector(sapply(test3, function(x) unlist(strsplit(x, "_"))[1]))
test4
test4[2] <- 'MSKPCa18'
test4[3] <- 'MSKPCa19'
test4[c(31,32,33)] <- c('MSKPCa20','MSKPCa22','MSKPCa24')
test4[4] <- 'C4-2'
test4[27:30] <- c("WCM1078","WCM1262","WCM154", "WCM155")
test4
names(organoid_peaks_all) <- test4

sampleReordered <- organoid_peaks_all[c("MSKPCa1", "MSKPCa2", "MSKPCa3", "MSKPCa8", "MSKPCa9", "MSKPCa10", "MSKPCa11", "MSKPCa12", "MSKPCa13", "MSKPCa14", "MSKPCa15", "MSKPCa16", "MSKPCa17", "MSKPCa18", "MSKPCa19", "MSKPCa20", "MSKPCa22", "MSKPCa24", "WCM154", "WCM155", "WCM1078", "WCM1262", "LNCaP", "VCaP", "22Rv1", "C4-2", "PC3", "DU145", "H660", "MSKEF1", "PARCB1", "PARCB3", "PARCB6", "PARCB8")]


peakAnnoList <- lapply(sampleReordered, annotatePeak, 
                       tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Hs.eg.db", TxDb=txdb)
saveRDS(peakAnnoList, 'figure1d_peakAnnoList_noPDX.RDS')

pdf("figure1/figure1d_noPDX.pdf", width=6, height=6, useDingbats=F)
plotAnnoBar(peakAnnoList,
            title = "Feature Distribution")
dev.off()


## feature distribution idr peaks w PDX
peaks <- dir("../EX_0201_prostateSubtypePredict/IDR_peaks/", pattern = "*.narrowPeak$", 
             full.names = TRUE)
organoid_peaks_all <- as.list(peaks)
organoid_peaks_all
test <- unlist(organoid_peaks_all)
test2 <- as.vector(sapply(test, function(x) unlist(strsplit(x, "IDR_peaks//"))[2]))
test3 <- as.vector(sapply(test2, function(x) unlist(strsplit(x, ".idr"))[1]))
test4 <- as.vector(sapply(test3, function(x) unlist(strsplit(x, "_"))[1]))
test4
test4[2] <- 'MSKPCa18'
test4[3] <- 'MSKPCa19'
test4[c(40,41,42)] <- c('MSKPCa20','MSKPCa22','MSKPCa24')
test4[4] <- 'C4-2'
test4[36:39] <- c("WCM1078","WCM1262","WCM154", "WCM155")
test4
names(organoid_peaks_all) <- test4

sampleReordered <- organoid_peaks_all[c("MSKPCa1", "MSKPCa2", "MSKPCa3", "MSKPCa8", "MSKPCa9", "MSKPCa10", "MSKPCa11", "MSKPCa12", "MSKPCa13", "MSKPCa14", "MSKPCa15", "MSKPCa16", "MSKPCa17", "MSKPCa18", "MSKPCa19", "MSKPCa20", "MSKPCa22", "MSKPCa24", "WCM154", "WCM155", "WCM1078", "WCM1262", "LNCaP", "VCaP", "22Rv1", "C4-2", "PC3", "DU145", "H660","PDX01a","PDX07b","PDX09a","PDX10a","PDX16a","PDX16b","PDX27a","PDX34a","PDX43a", "MSKEF1", "PARCB1", "PARCB3", "PARCB6", "PARCB8")]

peakAnnoList <- lapply(sampleReordered, annotatePeak, 
                       tssRegion=c(-3000, 3000), verbose=FALSE, annoDb="org.Hs.eg.db", TxDb=txdb)
saveRDS(peakAnnoList, 'figure1d_peakAnnoList_wPDX.RDS')

pdf("figure1/figure1d_wPDX.pdf", width=6, height=6, useDingbats=F)
plotAnnoBar(peakAnnoList,
            title = "Feature Distribution")
dev.off()
