library(DESeq2)
library(limma)
library(ConsensusClusterPlus)
library(pheatmap)
library(CancerSubtypes)
library(edgeR)

########################################################################################################################
# getting the top variable peaks using counts with replica
countData_wRep <- read.csv('featureCountMatrix_idr_all.csv', header = T, sep = ',', row.names = 1)
metaData_wRep <- read.csv('orgn_sampleList_wRep.txt', header = F, sep = "\t", row.names = 1)
consensusToCount <- readRDS('consensusToCount_idr.RDS')

colnames(countData_wRep)[c(1,2)]
colnames(countData_wRep)[1] <- '22Rv1_rep1'
colnames(countData_wRep)[2] <- '22Rv1_rep2'

dds_wRep <- DESeqDataSetFromMatrix(countData_wRep, metaData_wRep, ~ V2, rowRanges = consensusToCount, tidy = F)
dds_wRep <- DESeq(dds_wRep)

# log
normCounts_wRep <- counts(dds_wRep, normalized =TRUE)
rlog_wRep <- log2(normCounts_wRep+1)

# top % variable peaks function
top_var_gene_var <- function(input_matrix, threshold) {
  select <- head(order(rowVars(input_matrix),decreasing=TRUE),threshold)
  topVarGenes <- input_matrix[select, ]
  return(topVarGenes)
}

peak_topvar <- top_var_gene_var(rlog_wRep, as.integer(dim(rlog_wRep)[1]/100))

sample_cor <- cor(peak_topvar, method = "pearson")
annotation = data.frame(factor(metaData_wRep$V6))
rownames(annotation) = rownames(metaData_wRep)
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE, annotation = annotation)


########################################################################################################################
# do clustering using counts without replica
countData <- read.csv('featureCountMatrix_idr_ave.csv', header = T, sep = ',', row.names = 1)
metaData <- read.csv('../EX_0201_prostateSubtypePredict/orgn_sampleList_wPDX.txt', header = F, sep = "\t", row.names = 1)

colnames(countData)[1]
colnames(countData)[1] <- '22Rv1'

dds <- DESeqDataSetFromMatrix(countData=countData, colData = metaData, rowRanges = consensusToCount,
                              tidy = F, design=~1)

dds <- DESeq(dds)

# log
normCounts <- counts(dds, normalized =TRUE)
rlog <- log2(normCounts+1)

rlog <- rlog[,rownames(metaData)]

rlog_ave_topVar <- rlog[row.names(peak_topvar),]
rlog_ave_topVar <- sweep(rlog_ave_topVar,1, apply(rlog_ave_topVar,1,median,na.rm=T))
dim(rlog_ave_topVar)


library(pheatmap)
annotation = data.frame(factor(metaData$V6))
rownames(annotation) = rownames(metaData)
sample_cor <- cor(rlog_ave_topVar, method = "pearson")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE, annotation = annotation)
sample_cor <- cor(rlog_ave_topVar, method = "spearman")
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE, annotation = annotation)
# dev.off()


rlog_ave_topVar <- rlog[row.names(peak_topvar),]
rlog_ave_topVar <- sweep(rlog_ave_topVar,1, apply(rlog_ave_topVar,1,median,na.rm=T))

# Use top 1% variable peaks from with repeats (log2(norm+1)), and then center by median for clustering
# euclidean distance, hc, pItem anything > 0.98
results = ConsensusClusterPlus(as.matrix(rlog_ave_topVar),maxK=6,reps=1000,pItem=0.98,
                               pFeature=1,innerLinkage="ward.D2", finalLinkage="ward.D2",
                               title="consensusCluster_ave", distance='euclidean', writeTable=T,
                               clusterAlg="hc",plot="pdf",seed=42)

sample_cor <- cor(rlog_ave_topVar[,c(1:40)], method = "pearson")
dt = as.matrix(1-sample_cor)
pheatmap(sample_cor, cluster_rows=TRUE, show_rownames=TRUE, annotation = annotation)

calcRes <- calcICL(results, title = "consensusCluster_ave", plot = 'pdf')
write.csv(calcRes$clusterConsensus, 'consensusCluster_ave/clusterConsensus.csv', quote=F, row.names = F)
write.csv(calcRes$itemConsensus, 'consensusCluster_ave/itemConsensus.csv', quote=F, row.names = F)

ccLabel <- data.frame(factor(results[[4]]$consensusClass))
colnames(ccLabel) <- 'ccPredict'
ccMatrix <- results[[4]]$consensusMatrix
rownames(ccMatrix) <- names(results[[4]]$consensusClass)
colnames(ccMatrix) <- names(results[[4]]$consensusClass)

#pdf("consensusCluster_ave/ccHeatmap.pdf", width = 6.5, height = 5)
pheatmap(ccMatrix, cluster_rows=results[[4]]$consensusTree, cluster_cols=results[[4]]$consensusTree,
         show_rownames=TRUE, annotation_row = ccLabel, annotation_col = annotation, fontsize=6)
dev.off()

