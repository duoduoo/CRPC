library(DESeq2)
library(limma)
library(ConsensusClusterPlus)
library(pheatmap)
library(CancerSubtypes)
library(edgeR)

# top % variable peaks function
top_var_gene_var <- function(input_matrix, threshold) {
  select <- head(order(rowVars(input_matrix),decreasing=TRUE),threshold)
  topVarGenes <- input_matrix[select, ]
  return(topVarGenes)
}


########################################################################################################################
# do clustering using counts without replica
countData <- read.csv('featureCountMatrix_idr_ave.csv', header = T, sep = ',', row.names = 1)
metaData <- read.csv('orgn_sampleList_wPDX.txt', header = F, sep = "\t", row.names = 1)
consensusToCount <- readRDS('consensusToCount_idr.RDS')

colnames(countData)[1]
colnames(countData)[1] <- '22Rv1'


dds <- DESeqDataSetFromMatrix(countData=countData, colData = metaData, rowRanges = consensusToCount,
                              tidy = F, design=~1)

dds <- DESeq(dds)

# log
normCounts <- counts(dds, normalized =TRUE)
rlog <- log2(normCounts+1)

rlog <- rlog[,rownames(metaData)]

peak_topvar <- top_var_gene_var(rlog, as.integer(dim(rlog)[1]/100))

rlog_ave_topVar <- rlog[row.names(peak_topvar),]
rlog_ave_topVar <- sweep(rlog_ave_topVar,1, apply(rlog_ave_topVar,1,median,na.rm=T))

# Use top 1% variable peaks, (log2(norm+1)), and then center by median for clustering
results = ConsensusClusterPlus(as.matrix(rlog_ave_topVar),maxK=6,reps=1000,pItem=0.9,
                               pFeature=1,innerLinkage="ward.D", finalLinkage="ward.D2",
                               title="consensusCluster_ave", distance='pearson', writeTable=T,
                               clusterAlg="hc",plot="pdf",seed=23445.61)

ccLabel <- data.frame(factor(results[[4]]$consensusClass))
colnames(ccLabel) <- 'ccPredict'
ccMatrix <- results[[4]]$consensusMatrix
rownames(ccMatrix) <- names(results[[4]]$consensusClass)
colnames(ccMatrix) <- names(results[[4]]$consensusClass)

#pdf("consensusCluster_ave/ccHeatmap.pdf", width = 6.5, height = 5)
pheatmap(ccMatrix, cluster_rows=results[[4]]$consensusTree, cluster_cols=results[[4]]$consensusTree,
         show_rownames=TRUE, annotation_row = ccLabel, annotation_col = annotation, fontsize=6)
dev.off()
        
         
calcRes <- calcICL(results, title = "consensusCluster_ave", plot = 'pdf')
write.csv(calcRes$clusterConsensus, 'consensusCluster_ave/clusterConsensus.csv', quote=F, row.names = F)
write.csv(calcRes$itemConsensus, 'consensusCluster_ave/itemConsensus.csv', quote=F, row.names = F)




