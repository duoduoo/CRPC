
library(DESeq2)
library(limma)
library(ConsensusClusterPlus)
library(pheatmap)


countData <- read.csv('mapToHumanOnly/totalDf.csv', header = T, sep = ',', row.names = 1)
metaData <- read.csv('../EX_0201_prostateSubtypePredict/orgn_sampleList_wPDX.txt', header = F, sep = "\t", row.names = 1)


colnames(countData)[23]
colnames(countData)[23] <- '22Rv1'
countData <- countData[,row.names(metaData)]

countData <- countData[,1:40]
metaData <- metaData[1:40,]


dds <- DESeqDataSetFromMatrix(countData=countData, colData = metaData,
                              tidy = F, design=~ 1)


dds <- dds[rowSums(counts(dds) >0) > 10, ]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

sizeFacterOut <- matrix(dds$sizeFactor)
rownames(sizeFacterOut) <- names(dds$sizeFactor)
write.csv(sizeFacterOut, 'dds_sizeFactor.csv', quote=F)

# rlog
normCounts <- counts(dds, normalized =TRUE)
rlog <- log2(normCounts+1)
# vst
rlog <- vst(dds, blind=FALSE)
rlog<- assay(rlog)

# top 10% variable peaks function
top_var_gene_var <- function(input_matrix, threshold) {
  select <- head(order(rowVars(input_matrix),decreasing=TRUE),threshold)
  topVarGenes <- input_matrix[select, ]
  return(topVarGenes)
}


gene_topvar <- top_var_gene_var(rlog, 1000)



#########################################################################################################
# Plot RNA-seq heatmap figure 2G

# # ##### Get gene id - start ###############
library(biomaRt)
listMarts(host = 'grch37.ensembl.org')
ensembl <- useMart(host='grch37.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl',ensemblRedirect = F)

BiocManager::install("GSVA")
BiocManager::install("circlize")
library(GSVA)
library(circlize)
library(ComplexHeatmap)

head(rlog)

AR_sigName <- c("PTGER4","FKBP5","KLK2","CENPN","MAF","ACSL3","HERC3","ZBTB10","EAF2","ABCC4","C1orf116","PMEPA1","MED28","MPHOSPH9","TMPRSS2","KLK3","NKX3-1","NNMT","ADAM7","ELL2")
AR_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=AR_sigName, mart= ensembl, useCache = FALSE)

NEPC_sigName <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9","BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1","PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11","AURKA","DNMT1","EZH2","MYCN")
NEPC_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=NEPC_sigName, mart= ensembl, useCache = FALSE)

WNT_sigName <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2","FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2","HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1","MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4","NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53","WNT1","WNT5B","WNT6")
WNT_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=WNT_sigName, mart= ensembl, useCache = FALSE)

SCL_sigName <- c("ABI3BP","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","AEBP1","AGPAT4","AHI1","AKAP2","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARNTL","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","BMP1","BMP7","BNC1","BTBD11","BVES","C10orf54","C12orf24","C13orf33","C14orf149","C14orf37","C17orf81","C19orf12","C1QTNF4","C21orf7","C2orf40","C3orf54","C3orf58","C3orf64","C4orf49","C5orf25","C6orf145","C6orf228","C8orf84","C9orf3","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CBLB","CCDC3","CCDC85B","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CHSY3","CLIP3","CLMP","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTGF","CTNNAL1","CXCL14","CYGB","D4S234E","DCBLD2","DCHS1","DCUN1D3","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","EDARADD","EDNRB","EEPD1","EFCAB1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ENC1","ENPP2","EPAS1","EPDR1","EPHB1","ERF","ETS1","ETV5","EVC","EXT1","FABP5","FAM101B","FAM132A","FAM176A","FAM184A","FAM70B","FAS","FBLN1","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR124","GPR176","GPR3","GPR87","GPSM1","GRASP","GSN","GYLTL1B","GYPC","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KDELC1","KIAA0889","KLHDC5","KLHL21","KLHL29","KRT14","KRT16","KRT5","KRT75","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LBH","LCA5","LCAT","LEP","LEPRE1","LEPREL1","LGALS1","LGALS7","LGR6","LHFP","LIFR","LIMA1","LIMS2","LMOD1","LPHN1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MATN2","MBNL1","MCAM","MEF2C","MEG3","MEST","MFNG","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MRVI1","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NETO2","NGF","NGFR","NLGN2","NNAT","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NT5E","NTF3","NTRK2","NUDT10","NUDT11","NXN","ODZ3","OSBPL6","OSR1","OXTR","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PELO","PGF","PHLDA3","PHLDB1","PKD1","PKD2","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLS3","PLXNA2","PODN","POPDC2","POSTN","POU3F1","PPAP2A","PPAP2B","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPLA","PTPRE","PTPRT","PVRL3","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFX2","RGNEF","RHOJ","RND3","RNF165","RUSC2","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SDPR","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMTN","SNAI2","SNCA","SNTB2","SOBP","SORBS1","SORCS1","SOX11","SPARC","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TBX2","TCF4","TCF7L1","TCOF1","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178","TMEM201","TMEM204","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTL","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423")
SCL_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=SCL_sigName, mart= ensembl, useCache = FALSE)

gsvaScore <- gsva(rlog, list(AR_sigID$ensembl_gene_id_version, WNT_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version, SCL_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=15, max.sz=500, method="gsva")
gsvaScore
metaData[,c("AR_score","WNT_score","NEPC_score","SCL_score")] <- t(gsvaScore)[rownames(metaData),]

AR_list <- c("KLK2", "FKBP5", "TMPRSS2", "KLK3", "NKX3-1", "PMEPA1", "PPAP2A", "PART1", "ALDH1A3") #10:10
NE_list <- c("CHGA", "CHGB", "SCG3", "CHRNB2", "ELAVL4", "SYP", "ENO2", "SCN3A") #10:10
basal_list <- c("TP63", "TRIM29", "ITGB4", "KRT5", "KRT14") #5:5
luminal_list <- c("EPCAM", "KRT8", "KRT18", "HPN", "DPP4") #7:7
stem_list <- c("CD44","TACSTD2" ,"ATXN1") #8:8
FGFR_targets_list <- c("DUSP6", "SPRY4", "ETV4", "GDF15", "EGR1", "ETV5", "C3orf52", "UPP1", "DUSP5", "TRIB3", "MYEOV", "MAFF", "MESDC1", "ARPC5L", "MAP2K3", "YRDC", "IER3", "PHLDA1", "PMAIP1", "SPRED1", "LRP8") #21:21
YAP_TAZ_list <- c("YAP1", "WWTR1","CYR61", "CTGF", "AMOTL2", "ANKRD1", "F3","CRIM1", "AXL", "DOCK5", "MYOF", "AJUBA") #10:10
Wnt_list <- c("AXIN2", "RNF43", "TCF7L2", "CTNNB1") #4:4
EMT_list <- c("CDH1", "CDH2", "VIM", "ITGB4", "SNAI1", "SNAI2", "ZEB1", "TWIST1", "FN1", "ACTA2")
FGF_list <- c("FGF1", "FGF2", "FGF3", "FGF4", "FGF5", "FGF6", "FGF7", "FGF8", "FGF9", "FGF10", "FGF17", "FGF18", "FGF19", "FGF20", "FGF21", "FGF22", "FGF23", "FGFBP1")
FGFR_list <- c("FGFR1", "FGFR2", "FGFR3", "FGFR4")

# plot heatmap with selective gene expression
gene_select <- c("AR",AR_list, Wnt_list, NE_list, stem_list, basal_list, luminal_list)
gene_select_id <- getBM(filters= "external_gene_name",
                        attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                        values=gene_select, mart= ensembl, useCache = FALSE)
gene_select_id <- subset(gene_select_id, !startsWith(gene_select_id$chromosome_name, 'HG'))

rlog_heatmap <- rlog[match(gene_select_id$ensembl_gene_id_version, rownames(rlog)), rownames(metaData)[1:40]]

rownames(rlog_heatmap) <- gene_select_id$external_gene_name

colnames(metaData)[5] <- "ATAC-seq_group"
ann_colors = list(`ATAC-seq_group` = c(AR_dependent = "#F8766D", WNT = "#7CAE00", NEPC = "#00BFC4", SCL= "#C77CFF"))


ar_sample <- rownames(metaData[(metaData$`ATAC-seq_group` == 'AR_dependent'),])
ar_sample
ar_sample <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX_X0009aS1p3","PDX_X0010aS2p2","PDX_X0016aS2p1","PDX_X0016bS2p1","BM111","ST278","PDX_X0034aS1p4")
wnt_sample <- rownames(metaData[(metaData$`ATAC-seq_group` == 'WNT'),])
wnt_sample
wnt_sample <- c("MSKPCa1","MSKPCa16","PM1078","PM1262")
nepc_sample <- rownames(metaData[(metaData$`ATAC-seq_group` == 'NEPC'),])
nepc_sample
nepc_sample <- c("MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","PM154","H660","ST280","PRNEX01aS1")
scl_sample <- rownames(metaData[(metaData$`ATAC-seq_group` == 'SCL'),])
scl_sample
scl_sample <- c("MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","PM155","DU145","PC3","BM110","ST262")


sample_order <- c(ar_sample, wnt_sample, nepc_sample, scl_sample)
length(sample_order)

# rlog_heatmap[which(rownames(rlog_heatmap)=="KRT5"),][28] <-2.571942
mat_scaled <-  t(apply(rlog_heatmap,1,scale)) 
colnames(mat_scaled) <- colnames(rlog_heatmap)
rownames(mat_scaled) <- rownames(rlog_heatmap)
mat_scaled <- as.matrix(mat_scaled)
# column_labels_renamed <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX09a","PDX10a","PDX16a","PDX16b","MSKPCa19","MSKPCa22","PDX34a","MSKPCa1","MSKPCa16","WCM1078","WCM1262","PDX27a","MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","WCM154","H660","MSKPCa24","PDX01a","PDX07b","PDX43a","MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","WCM155","DU145","PC3","MSKPCa18","MSKPCa20")
sample_order
column_labels_renamed <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX09a","PDX10a","PDX16a","PDX16b","MSKPCa19","MSKPCa22","PDX34a","MSKPCa1","MSKPCa16","WCM1078","WCM1262","MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","WCM154","H660","MSKPCa24","PDX01a","MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","WCM155","DU145","PC3","MSKPCa18","MSKPCa20")
column_labels_renamed
length(column_labels_renamed)

col_AR <- colorRamp2(c(quantile(metaData$AR_score, 0), quantile(metaData$AR_score, 1)), c("white","#F8766D"))
col_WNT <- colorRamp2(c(quantile(metaData$WNT_score, 0), quantile(metaData$WNT_score, 1)), c("white","#7CAE00"))
col_NE <- colorRamp2(c(quantile(metaData$NEPC_score, 0), quantile(metaData$NEPC_score, 1)), c("white","#00BFC4"))
col_CSC <- colorRamp2(c(quantile(metaData$SCL_score, 0), quantile(metaData$SCL_score, 1)), c("white","#C77CFF"))

ha_top <- HeatmapAnnotation(df=metaData[sample_order, c("ATAC-seq_group", "AR_score","WNT_score", "NEPC_score","SCL_score")],
                            col = list(`ATAC-seq_group`=c("AR_dependent" = "#F8766D", "WNT" = "#7CAE00", "NEPC" = "#00BFC4", "SCL"= "#C77CFF"),
                                       AR_score=col_AR,
                                       WNT_score=col_WNT,
                                       NEPC_score=col_NE,
                                       SCL_score=col_CSC), 
                            show_annotation_name = TRUE)

ht_list <- Heatmap(mat_scaled[gene_select,sample_order], 
                   top_annotation = ha_top,  
                   cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = T,
                   column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 10),
                   rect_gp = gpar(col = "black", lwd = 0.3),
                   column_labels = column_labels_renamed,
                   column_title = "Selective gene expression in organoids and cell lines", 
                   column_title_gp = gpar(fontsize = 14, fontface = "bold"))

draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "left")


dev.off()


#########################################################################################################
# Plot RNA-seq heatmap figure S3B

sigGenes <- read.csv('geneSignature_noPDX_noBatch_cutoff_05.csv')
ar_gene <- sigGenes[(sigGenes$class == 'AR_dependent'),'probe']
wnt_gene <- sigGenes[(sigGenes$class == 'WNT'),'probe']
ne_gene <- sigGenes[(sigGenes$class == ''),'probe']
scl_gene <- sigGenes[(sigGenes$class == 'AR_dependent'),'probe']


ar_sample <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX_X0009aS1p3","PDX_X0010aS2p2","PDX_X0016aS2p1","PDX_X0016bS2p1","BM111","ST278","PDX_X0034aS1p4")
wnt_sample <- c("MSKPCa1","MSKPCa16","PM1078","PM1262")
nepc_sample <- c("MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","PM154","H660","ST280","PRNEX01aS1")
scl_sample <- c("MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","PM155","DU145","PC3","BM110","ST262")
sample_order <- c(ar_sample, wnt_sample, nepc_sample, scl_sample)


rlog_reorder <- rlog[match(sigGenes$probe, rownames(rlog)), match(sample_order, colnames(rlog))]

pdf(paste0("figures/", "unique_sign_heatmap-", Sys.Date(), ".pdf"), width=5, height=6, useDingbats=F)
pheatmap(rlog_reorder, cluster_rows=F, cluster_cols = F, fontsize_row = 0.4, fontsize_col = 9, scale = "row",
         annotation_col = metaData[, c(5,9)],
         main="upregulated genes in each subtype")
dev.off()


#########################################################################################################
# Plot RNA-seq heatmap figure S4B

gene_select <- FGFR_targets_list
gene_select_id <- getBM(filters= "external_gene_name",
                        attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                        values=gene_select, mart= ensembl, useCache = FALSE)
gene_select_id <- subset(gene_select_id, !startsWith(gene_select_id$chromosome_name, 'H'))

gsvaScore <- gsva(as.matrix(rlog), list(gene_select_id$ensembl_gene_id_version), mx.diff=F, min.sz=15, max.sz=500, method="gsva")
gsvaScore

metaData <- metaData[1:40,]

metaData[,c("FGF_score")] <- t(gsvaScore)[rownames(metaData),]
saveRDS(metaData, 'rnaSeq_metaData_drop3pdx.RDS')
metaData <- readRDS('rnaSeq_metaData_drop3pdx.RDS')

ar_sample <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX_X0009aS1p3","PDX_X0010aS2p2","PDX_X0016aS2p1","PDX_X0016bS2p1","BM111","ST278","PDX_X0034aS1p4")
wnt_sample <- c("MSKPCa1","MSKPCa16","PM1078","PM1262")
nepc_sample <- c("MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","PM154","H660","ST280","PRNEX01aS1")
scl_sample <- c("MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","PM155","DU145","PC3","BM110","ST262")

sample_order <- c(ar_sample, wnt_sample, nepc_sample, scl_sample)
length(sample_order)

FGFR_targets_list_id <- getBM(filters= "external_gene_name",
                              attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                              values=FGFR_targets_list, mart= ensembl, useCache = FALSE)
FGFR_targets_list_id <- subset(FGFR_targets_list_id, !startsWith(FGFR_targets_list_id$chromosome_name, 'H'))

FGF_list_id <- getBM(filters= "external_gene_name",
                              attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                              values=FGF_list, mart= ensembl, useCache = FALSE)
FGF_list_id <- subset(FGF_list_id, !startsWith(FGF_list_id$chromosome_name, 'H'))

# FGF genes
rlog_heatmap <- rlog[match(FGF_list_id$ensembl_gene_id_version, rownames(rlog)), sample_order]
rownames(rlog_heatmap) <- FGF_list_id$external_gene_name

# FGF target genes 
rlog_heatmap <- rlog[match(FGFR_targets_list_id$ensembl_gene_id_version, rownames(rlog)), sample_order]
rownames(rlog_heatmap) <- FGFR_targets_list_id$external_gene_name

dim(rlog_heatmap)
colnames(metaData)[5] <- "ATAC-seq_group"
ann_colors = list(`ATAC-seq_group` = c(AR_dependent = "#F8766D", WNT = "#7CAE00", NEPC = "#00BFC4", SCL= "#C77CFF"))


mat_scaled <-  t(apply(rlog_heatmap,1,scale)) 
colnames(mat_scaled) <- colnames(rlog_heatmap)
rownames(mat_scaled) <- rownames(rlog_heatmap)
mat_scaled <- as.matrix(mat_scaled)

sample_order
column_labels_renamed <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX09a","PDX10a","PDX16a","PDX16b","MSKPCa19","MSKPCa22","PDX34a","MSKPCa1","MSKPCa16","WCM1078","WCM1262","MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","WCM154","H660","MSKPCa24","PDX01a","MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","WCM155","DU145","PC3","MSKPCa18","MSKPCa20")
column_labels_renamed
length(column_labels_renamed)

pheatmap(mat_scaled[FGF_list,sample_order], cluster_rows=F, cluster_cols = F, fontsize_row = 8, fontsize_col = 6, scale = "row", show_rownames = T,
         annotation_col = metaData[sample_order, c("ATAC-seq_group", "FGF_score","AR_score","NEPC_score")],labels_col = column_labels_renamed,
         annotation_colors = ann_colors,main="Expression of FGF ligands and receptors")

pheatmap(mat_scaled[FGFR_targets_list,sample_order], cluster_rows=F, cluster_cols = F, fontsize_row = 8, fontsize_col = 6, scale = "row", show_rownames = T,
         annotation_col = metaData[sample_order, c("ATAC-seq_group", "FGF_score","AR_score","NEPC_score")],labels_col = column_labels_renamed,
         annotation_colors = ann_colors,main="Expression of FGF targets")


dev.off()

#########################################################################################################
# Plot RNA-seq heatmap figure S14A, AP-1 component genes

gene_select <- c('FOS','FOSB','FOSL1', 'FOSL2','JUN','JUNB','JUND')
gene_select_id <- getBM(filters= "external_gene_name",
                        attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                        values=gene_select, mart= ensembl, useCache = FALSE)
gene_select_id <- subset(gene_select_id, !startsWith(gene_select_id$chromosome_name, 'H'))

metaData <- readRDS('rnaSeq_metaData_drop3pdx.RDS')

ar_sample <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX_X0009aS1p3","PDX_X0010aS2p2","PDX_X0016aS2p1","PDX_X0016bS2p1","BM111","ST278","PDX_X0034aS1p4")
wnt_sample <- c("MSKPCa1","MSKPCa16","PM1078","PM1262")
nepc_sample <- c("MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","PM154","H660","ST280","PRNEX01aS1")
scl_sample <- c("MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","PM155","DU145","PC3","BM110","ST262")

sample_order <- c(ar_sample, wnt_sample, nepc_sample, scl_sample)
length(sample_order)

# AP-1 component genes
rlog_heatmap <- rlog[match(gene_select_id$ensembl_gene_id_version, rownames(rlog)), sample_order]
rownames(rlog_heatmap) <- gene_select_id$external_gene_name

mapply(sum, rlog_heatmap)

dim(rlog)
gsvaScore <- gsva(as.matrix(rlog), list(gene_select_id$ensembl_gene_id_version), mx.diff=F, min.sz=5, max.sz=500, method="gsva")
gsvaScore
View(metaData)
metaData[,c("AP1_score")] <- t(gsvaScore)[rownames(metaData),]
metaData[,c("AP1_score")] <-mapply(sum, rlog_heatmap)


dim(rlog_heatmap)
colnames(metaData)[5] <- "ATAC-seq_group"
ann_colors = list(`ATAC-seq_group` = c(AR_dependent = "#F8766D", WNT = "#7CAE00", NEPC = "#00BFC4", SCL= "#C77CFF"))

# rlog_heatmap[which(rownames(rlog_heatmap)=="KRT5"),][28] <-2.571942
mat_scaled <-  t(apply(rlog_heatmap,1,scale))
colnames(mat_scaled) <- colnames(rlog_heatmap)
rownames(mat_scaled) <- rownames(rlog_heatmap)
mat_scaled <- as.matrix(mat_scaled)

column_labels_renamed <- c("MSKPCa2","22Rv1","C4.2","LNCaP","VCaP","PDX09a","PDX10a","PDX16a","PDX16b","MSKPCa19","MSKPCa22","PDX34a","MSKPCa1","MSKPCa16","WCM1078","WCM1262","MSKEF1","PARCB1","PARCB3","PARCB6","PARCB8","MSKPCa10","MSKPCa14","WCM154","H660","MSKPCa24","PDX01a","MSKPCa11","MSKPCa12","MSKPCa13","MSKPCa15","MSKPCa17","MSKPCa3","MSKPCa8","MSKPCa9","WCM155","DU145","PC3","MSKPCa18","MSKPCa20")
column_labels_renamed
length(column_labels_renamed)

pheatmap(mat_scaled[gene_select,sample_order], cluster_rows=F, cluster_cols = F, fontsize_row = 10, fontsize_col = 10, scale = "row", show_rownames = T,
         annotation_col = metaData[sample_order, c("ATAC-seq_group", "AP1_score")],labels_col = column_labels_renamed,
         annotation_colors = ann_colors, main="Expression of AP-1 components",
         color=colorRampPalette(c("blue","white", "red"))(30))


dev.off()

