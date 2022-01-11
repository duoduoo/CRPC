setwd("~/Documents/KhuranaLab/EX_0209_crpcFigures")

BiocManager::install("survival")
BiocManager::install("survminer")

library(Biobase)
library(CMScaller)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(survival)
library(survminer)
library(dplyr)
library(ranger)
library(ggplot2)
library(ggfortify)
#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(GSVA)
library(ComplexHeatmap)
library(survival)
library(survminer)


####################################################################################################################################
# Table_S8 signature genes (hg19)
tmp_org <- read.csv('../EX_0208_prostateRnaCluster/geneSignature_noPDX_noBatch.csv', header = T)
tmp_org <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/rankBy_baseMean_top93.csv', header = T)

# # ##### Get gene id - start ###############
listMarts(host = 'grch37.ensembl.org')
ensembl <- useMart(host='grch37.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')

gene_name <- getBM(filters= "ensembl_gene_id_version",
                   attributes= c("ensembl_gene_id","ensembl_gene_id_version","chromosome_name", "external_gene_name", "gene_biotype"),
                   values=tmp_org$probe, mart= ensembl, useCache = FALSE)
row.names(gene_name) <- gene_name$ensembl_gene_id_version

head(gene_name)
tmp_org$gene_name <- gene_name[match(tmp_org$probe, rownames(gene_name)),'external_gene_name']
write.csv(tmp_org, 'tables/table_S8.csv', quote=F, row.names = F)

####################################################################################################################################
# figure 4B. IPM/SU2C patient--NTP plot

# listMarts(host = 'uswest.ensembl.org')
ensembl_hg38 <- useMart(host='uswest.ensembl.org',
                   biomart='ENSEMBL_MART_ENSEMBL',
                   dataset='hsapiens_gene_ensembl')

tmp_org <- read.csv('../EX_0208_prostateRnaCluster/geneSignature_noPDX_noBatch.csv', header = T)
tmp_org <- read.csv('../EX_0208_prostateRnaCluster/geneSignature_noPDX_noBatch_hg38_su2c.csv', header = T)
tmp_org <- tmp_org[c(which(tmp_org$class == 'AR_dependent'),which(tmp_org$class == 'WNT'),which(tmp_org$class == 'NEPC'),which(tmp_org$class == 'SCL')),]

####################################################################################################################################
# # not used
# patient_combat_data <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/IPM_patient_RSEM.csv', row.names = 2)
# patient_combat_data <- patient_combat_data[,-c(1)]
# patient_combat_data <- read.csv('../EX_0208_prostateRnaCluster/su2c_rnaseq_geneNameRevised_noPDX_noBatch.csv', header = T, row.names = 1)
# patient_combat_data <- read.csv('../EX_0208_prostateRnaCluster/su2c_tpm_geneID_fill0.csv', header = T, row.names = 1)
####################################################################################################################################

#patient_combat_data <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/IPM_patient_RSEM_combine3MSK.csv', row.names = 1)
patient_combat_data <- read.csv('../EX_0208_prostateRnaCluster/su2c_rnaseq_geneNameRevised_noPDX_noBatch_266Sample_new.csv', header = T, row.names = 1)

#patient_combat_data_meta <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/IPM_patient_meta_Shaham_revised_combine3MSK.csv', sep=',',row.names = 1,fileEncoding = "latin1")
patient_combat_data_meta <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/pnas.1902651116.sd01_revised.csv', row.names = 2)

####patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_ipmPatients_noPDX_noBatch.csv', row.names = 1)
####patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_su2cPatients_noPDX_noBatch.csv', row.names = 1)

#patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/RES_rankBy_baseMean_top93_ipm.csv', row.names = 1)
patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/RES_rankBy_baseMean_top93_su2c.csv', row.names = 1)


patient_combat_data_result <- patient_combat_data_result[which(patient_combat_data_result$FDR < 0.5),]
patient_combat_data_result <- patient_combat_data_result[which(patient_combat_data_result$p.value < 0.05),]

ar_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'AR_dependent'),])
wnt_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'WNT'),])
nepc_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'NEPC'),])
scl_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'SCL'),])

patient_order <- c(ar_patient, wnt_patient, nepc_patient, scl_patient)

# re-order by group
patient_combat_data <- patient_combat_data[,match(patient_order, colnames(patient_combat_data))]
patient_combat_data_result <- patient_combat_data_result[match(patient_order, rownames(patient_combat_data_result)),]
patient_combat_data_meta <- patient_combat_data_meta[match(patient_order, rownames(patient_combat_data_meta)),]

# calculate AR and NE score for patients
rlog <- log2(patient_combat_data+1)


AR_sigName <- c("PTGER4","FKBP5","KLK2","CENPN","MAF","ACSL3","HERC3","ZBTB10","EAF2","ABCC4","C1orf116","PMEPA1","MED28","MPHOSPH9","TMPRSS2","KLK3","NKX3-1","NNMT","ADAM7","ELL2")
length(AR_sigName)
AR_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=AR_sigName, mart= ensembl, useCache = FALSE)
AR_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=AR_sigName, mart=ensembl_hg38, useCache = FALSE)
row.names(AR_sigID) <- AR_sigID$external_gene_name

NEPC_sigName <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9","BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1","PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11","AURKA","DNMT1","EZH2","MYCN")
length(NEPC_sigName)
NEPC_sigID <- getBM(filters= "external_gene_name",
                    attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                    values=NEPC_sigName, mart= ensembl, useCache = FALSE)
NEPC_sigID <- getBM(filters= "external_gene_name",
                    attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                    values=NEPC_sigName, mart= ensembl_hg38, useCache = FALSE)
row.names(NEPC_sigID) <- NEPC_sigID$external_gene_name


row.names(rlog)[which(row.names(rlog) %in% AR_sigID$external_gene_name)] <- AR_sigID[row.names(rlog)[which(row.names(rlog) %in% AR_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% NEPC_sigID$external_gene_name)] <- NEPC_sigID[row.names(rlog)[which(row.names(rlog) %in% NEPC_sigID$external_gene_name)],'ensembl_gene_id_version']

gsvaScore <- gsva(as.matrix(rlog), list(AR_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=15, max.sz=500, method="gsva")
# gsvaScore <- gsva(as.matrix(rlog), list(AR_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=10, max.sz=500, method="gsva")
gsvaScore
patient_combat_data_result[,c("AR_score","NEPC_score")] <- t(gsvaScore)[rownames(patient_combat_data_result),]


# select only the signature genes
mat <- rlog[match(tmp_org$probe, rownames(patient_combat_data)), match(patient_order, colnames(patient_combat_data))]
mat_scaled_2 <-  sweep(mat,1, apply(mat,1,mean,na.rm=T)) 
mat_scaled <-  t(apply(mat,1,scale)) 
colnames(mat_scaled) <- colnames(mat)
rownames(mat_scaled) <- rownames(mat)
mat_scaled <- as.data.frame(mat_scaled)


df_final <- patient_combat_data_result[, c(1,8,9)]
#df_final$pathology <- patient_combat_data_meta$Classification
df_final$pathology <- patient_combat_data_meta$Pathology.Classification
colnames(df_final) <- c("molecular_subtypes", "AR_score", "NE_score", "pathology")
df_final <- df_final[c("molecular_subtypes", "pathology", "NE_score", "AR_score")]

df_row = data.frame(tmp_org$class)
colnames(df_row)[1] = "subtype_signature"

col_AR <- colorRamp2(c(quantile(df_final$AR_score, 0), quantile(df_final$AR_score, 0.5), quantile(df_final$AR_score, 1)), c("green","white", "red"))
col_NE <- colorRamp2(c(quantile(df_final$NE_score, 0), quantile(df_final$NE_score, 0.5), quantile(df_final$NE_score, 1)), c("blue","white", "orange"))

#ha = HeatmapAnnotation(df = df_final, col = list(molecular_subtypes = c("AR_dependent" = "#F8766D", "WNT" = "#7CAE00", "NEPC"= "#00BFC4", "SCL"="#C77CFF", "class5_unknown"="grey"), pathology=c("CRPC-NEPC"="purple", "CRPC-Adeno"="pink", "CRPC-Squamous"="dark blue"), AR_score=col_AR, NE_score=col_NE))
ha = HeatmapAnnotation(df = df_final, col = list(molecular_subtypes = c("AR_dependent" = "#F8766D", "WNT" = "#7CAE00", "NEPC"= "#00BFC4", "SCL"="#C77CFF", "class5_unknown"="grey"), pathology = c("Adenocarcinoma" =  "pink", "Adenocarcinoma with NE features"= "white", "Small cell"="purple", "Not available"="grey", "Inadequate for diagnosis"="grey"), AR_score=col_AR, NE_score=col_NE))
ha_row = rowAnnotation(df = df_row, col = list(subtype_signature = c("AR_dependent" = "#F8766D", "WNT" = "#7CAE00", "NEPC"= "#00BFC4", "SCL"="#C77CFF", "class5_unknown"="grey")), show_annotation_name = F)

ht_list <- Heatmap(as.matrix(mat_scaled), 
                   top_annotation = ha,
                   cluster_rows = F,
                   cluster_columns = F,
                   show_row_names = F,
                   show_column_names = F,
                   show_heatmap_legend = T,
                   column_title = "Molecular classification of IPM patients", 
                   column_title_gp = gpar(fontsize = 14, fontface = "bold")) + ha_row

saveRDS(ht_list, paste0("figure4/", "IPM_NTP_results_scaled-", Sys.Date(), ".RDS"))
pdf(paste0("figure4/", "IPM_NTP_results_scaled-", Sys.Date(), ".pdf"), width=6, height=6, useDingbats=F)

saveRDS(ht_list, paste0("figure4/", "SU2C_NTP_results_scaled-", Sys.Date(), ".RDS"))
# ht_list <- readRDS(paste0("figure4/", "SU2C_NTP_results_scaled-", Sys.Date(), ".RDS"))
pdf(paste0("figure4/", "SU2C_NTP_results_scaled-", Sys.Date(), ".pdf"), width=7, height=6, useDingbats=F)

draw(ht_list, heatmap_legend_side = "right", annotation_legend_side = "left")
dev.off()


##########################################################################################################################################################
##########################################################################################################################################################
# Figure S13

# IPM patient
patient_combat_data <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/IPM_patient_RSEM_combine3MSK.csv', row.names = 1)
# patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_ipmPatients_noPDX_noBatch.csv', row.names = 1)
patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/RES_rankBy_baseMean_top93_ipm.csv', row.names = 1)

# rlog <- patient_combat_data
rlog <- log2(patient_combat_data+1)
rlog <- as.data.frame(scale(rlog, center = T, scale = apply(rlog, 2, sd, na.rm = TRUE)))
#rlog <- as.data.frame(scale(rlog, center = T, scale = rep(sd(unlist(rlog)),length(colnames(rlog)))))

AR_sigName <- c("PTGER4","FKBP5","KLK2","CENPN","MAF","ACSL3","HERC3","ZBTB10","EAF2","ABCC4","C1orf116","PMEPA1","MED28","MPHOSPH9","TMPRSS2","KLK3","NKX3-1","NNMT","ADAM7","ELL2")
# only for s23d
AR_sigName <- c("GNMT", "PTGER4", "FKBP5", "KLK2", "CENPN", "MAF", "ACSL3", "HERC3", "ZBTB10", "EAF2", "ABCC4", "C1orf116", "PMEPA1", "MED28", "MPHOSPH9", "TMPRSS2", "KLK3", "NKX3-1", "NNMT", "ADAM7", "ELL2")
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
WNT_sigID <- subset(WNT_sigID, !startsWith(WNT_sigID$chromosome_name, 'H'))

SCL_sigName <- c("ABI3BP","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","AEBP1","AGPAT4","AHI1","AKAP2","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARNTL","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","BMP1","BMP7","BNC1","BTBD11","BVES","C10orf54","C12orf24","C13orf33","C14orf149","C14orf37","C17orf81","C19orf12","C1QTNF4","C21orf7","C2orf40","C3orf54","C3orf58","C3orf64","C4orf49","C5orf25","C6orf145","C6orf228","C8orf84","C9orf3","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CBLB","CCDC3","CCDC85B","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CHSY3","CLIP3","CLMP","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTGF","CTNNAL1","CXCL14","CYGB","D4S234E","DCBLD2","DCHS1","DCUN1D3","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","EDARADD","EDNRB","EEPD1","EFCAB1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ENC1","ENPP2","EPAS1","EPDR1","EPHB1","ERF","ETS1","ETV5","EVC","EXT1","FABP5","FAM101B","FAM132A","FAM176A","FAM184A","FAM70B","FAS","FBLN1","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR124","GPR176","GPR3","GPR87","GPSM1","GRASP","GSN","GYLTL1B","GYPC","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KDELC1","KIAA0889","KLHDC5","KLHL21","KLHL29","KRT14","KRT16","KRT5","KRT75","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LBH","LCA5","LCAT","LEP","LEPRE1","LEPREL1","LGALS1","LGALS7","LGR6","LHFP","LIFR","LIMA1","LIMS2","LMOD1","LPHN1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MATN2","MBNL1","MCAM","MEF2C","MEG3","MEST","MFNG","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MRVI1","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NETO2","NGF","NGFR","NLGN2","NNAT","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NT5E","NTF3","NTRK2","NUDT10","NUDT11","NXN","ODZ3","OSBPL6","OSR1","OXTR","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PELO","PGF","PHLDA3","PHLDB1","PKD1","PKD2","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLS3","PLXNA2","PODN","POPDC2","POSTN","POU3F1","PPAP2A","PPAP2B","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPLA","PTPRE","PTPRT","PVRL3","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFX2","RGNEF","RHOJ","RND3","RNF165","RUSC2","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SDPR","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMTN","SNAI2","SNCA","SNTB2","SOBP","SORBS1","SORCS1","SOX11","SPARC","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TBX2","TCF4","TCF7L1","TCOF1","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178","TMEM201","TMEM204","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTL","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423")
SCL_sigID <- getBM(filters= "external_gene_name",
                   attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                   values=SCL_sigName, mart= ensembl, useCache = FALSE)
SCL_sigID <- subset(SCL_sigID, !startsWith(SCL_sigID$chromosome_name, 'H'))

YAP_TAZ_list <- c("YAP1", "WWTR1","CYR61", "CTGF", "AMOTL2", "ANKRD1", "F3","CRIM1", "AXL", "DOCK5", "MYOF", "AJUBA") #10:10
basal_list <- c("TP63", "TRIM29", "ITGB4", "KRT5", "KRT14")

YAP_TAZ_list <- c('CYR61','CTGF','AMOTL2','ANKRD1','IGFBP3','F3','FJX1','NUAK2','LAST2','CRIM1','GADD45A','TGFB2','PTPN14','NTSE','FOXF2','AXL','DOCK5','ASAP1','RBMS3','MYOF','ANHGEF17','CCDC80','AJUBA')
YAP_TAZ_sigID <- getBM(filters= "external_gene_name",
                   attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                   values=YAP_TAZ_list, mart= ensembl, useCache = FALSE)
YAP_TAZ_sigID <- subset(YAP_TAZ_sigID, !startsWith(YAP_TAZ_sigID$chromosome_name, 'H'))


basal_list <- c('COL17A1','CSMD2','CDH13','MUM1L1','MMP3','IL33','GIMAP8','PDPN','VSNL1','BNC1','IGFBP7','DLK2','HMGA2','NOTCH4','THBS2','TAGLN','FHL1','ANXA8L2','COL4A6','KCNQ5','WNT7A','KCNMA1','NIPAL4','FLRT2','LTBP2','FOXI1','NGFR','SERPINB13','CNTNAP3B','FGFR3','ARHGAP25','AEBP1','FJX1','TNC','MSRB3','NRG1','SERPINF1','DLC1','IL1A','DKK3','ERG','SYNE1','JAG2','JAM3','MRC2','SPARC','C16orf74','FAT3','KIRREL','SH2D5','KRT6A','KRT34','ITGA6','TP63','KRT5','KRT14')
basal_sigID <- getBM(filters= "external_gene_name",
                       attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                       values=basal_list, mart= ensembl, useCache = FALSE)
basal_sigID <- subset(basal_sigID, !startsWith(basal_sigID$chromosome_name, 'H'))


gsvaScore <- gsva(as.matrix(rlog), list(AR_sigID$ensembl_gene_id_version, WNT_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version, SCL_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=5, max.sz=500, method="zscore")
gsvaScore
patient_combat_data_result[,c("AR_score","WNT_score","NEPC_score","SCL_score")] <- t(gsvaScore)[rownames(patient_combat_data_result),]

patient_combat_data_result[,c("AR_score","WNT_score","NEPC_score","SCL_score",'YAP_TAZ_score','basal_score')] <- cbind(colSums(rlog[AR_sigID$ensembl_gene_id_version,]),colSums(rlog[WNT_sigID$ensembl_gene_id_version,]),colSums(rlog[NEPC_sigID$ensembl_gene_id_version,]),colSums(rlog[SCL_sigID$ensembl_gene_id_version,]),colSums(rlog[YAP_TAZ_sigID$ensembl_gene_id_version,]),colSums(rlog[basal_sigID$ensembl_gene_id_version,]))

# write.csv(patient_combat_data_result, '../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_ipmPatients_noPDX_noBatch_wSignalScore.csv')
write.csv(patient_combat_data_result, '../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_ipmPatients_top93_wSignalScore.csv')


# SU2C patient
# patient_combat_data <- read.csv('../EX_0208_prostateRnaCluster/su2c_rnaseq_geneNameRevised_noPDX_noBatch_266Samples.csv', header = T, row.names = 1)
# patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_su2cPatients_noPDX_noBatch.csv', row.names = 1)

patient_combat_data <- read.csv('../EX_0208_prostateRnaCluster/su2c_rnaseq_geneNameRevised_noPDX_noBatch_266Sample_new.csv', header = T, row.names = 1)
patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/RES_rankBy_baseMean_top93_su2c.csv', row.names = 1)

rlog <- log2(patient_combat_data+1)
rlog <- as.data.frame(scale(rlog, center = T, scale = apply(rlog, 2, sd, na.rm = TRUE)))
# rlog <- as.data.frame(scale(rlog, center = T, scale = rep(sd(unlist(rlog)),length(colnames(rlog)))))

AR_sigName <- c("PTGER4","FKBP5","KLK2","CENPN","MAF","ACSL3","HERC3","ZBTB10","EAF2","ABCC4","C1orf116","PMEPA1","MED28","MPHOSPH9","TMPRSS2","KLK3","NKX3-1","NNMT","ADAM7","ELL2")
AR_sigID <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=AR_sigName, mart=ensembl_hg38, useCache = FALSE)
row.names(AR_sigID) <- AR_sigID$external_gene_name

# remove 'KIAA0408'
NEPC_sigName <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","SOGA3","TRIM9","BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1","PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11","AURKA","DNMT1","EZH2","MYCN")
NEPC_sigID <- getBM(filters= "external_gene_name",
                    attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                    values=NEPC_sigName, mart= ensembl_hg38, useCache = FALSE)
NEPC_sigID <- subset(NEPC_sigID, NEPC_sigID$gene_biotype == 'protein_coding')
row.names(NEPC_sigID) <- NEPC_sigID$external_gene_name

WNT_sigName <- c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2","FRAT1","FZD1","FZD8","GNAI1","HDAC11","HDAC2","HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1","MAML1","MYC","NCOR2","NCSTN","NKD1","NOTCH1","NOTCH4","NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53","WNT1","WNT5B","WNT6")
WNT_sigID <- getBM(filters= "external_gene_name",
                   attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                   values=WNT_sigName, mart= ensembl_hg38, useCache = FALSE)
WNT_sigID <- subset(WNT_sigID, !startsWith(WNT_sigID$chromosome_name, 'CH'))
row.names(WNT_sigID) <- WNT_sigID$external_gene_name

SCL_sigName <- c("ABI3BP","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","AEBP1","AGPAT4","AHI1","AKAP2","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARNTL","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","BMP1","BMP7","BNC1","BTBD11","BVES","C10orf54","C12orf24","C13orf33","C14orf149","C14orf37","C17orf81","C19orf12","C1QTNF4","C21orf7","C2orf40","C3orf54","C3orf58","C3orf64","C4orf49","C5orf25","C6orf145","C6orf228","C8orf84","C9orf3","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CBLB","CCDC3","CCDC85B","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CHSY3","CLIP3","CLMP","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTGF","CTNNAL1","CXCL14","CYGB","D4S234E","DCBLD2","DCHS1","DCUN1D3","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","EDARADD","EDNRB","EEPD1","EFCAB1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ENC1","ENPP2","EPAS1","EPDR1","EPHB1","ERF","ETS1","ETV5","EVC","EXT1","FABP5","FAM101B","FAM132A","FAM176A","FAM184A","FAM70B","FAS","FBLN1","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR124","GPR176","GPR3","GPR87","GPSM1","GRASP","GSN","GYLTL1B","GYPC","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KDELC1","KIAA0889","KLHDC5","KLHL21","KLHL29","KRT14","KRT16","KRT5","KRT75","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LBH","LCA5","LCAT","LEP","LEPRE1","LEPREL1","LGALS1","LGALS7","LGR6","LHFP","LIFR","LIMA1","LIMS2","LMOD1","LPHN1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MATN2","MBNL1","MCAM","MEF2C","MEG3","MEST","MFNG","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MRVI1","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NETO2","NGF","NGFR","NLGN2","NNAT","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NT5E","NTF3","NTRK2","NUDT10","NUDT11","NXN","ODZ3","OSBPL6","OSR1","OXTR","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PELO","PGF","PHLDA3","PHLDB1","PKD1","PKD2","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLS3","PLXNA2","PODN","POPDC2","POSTN","POU3F1","PPAP2A","PPAP2B","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPLA","PTPRE","PTPRT","PVRL3","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFX2","RGNEF","RHOJ","RND3","RNF165","RUSC2","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SDPR","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMTN","SNAI2","SNCA","SNTB2","SOBP","SORBS1","SORCS1","SOX11","SPARC","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TBX2","TCF4","TCF7L1","TCOF1","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178","TMEM201","TMEM204","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTL","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423")
length(SCL_sigName)
SCL_sigID <- getBM(filters= "external_gene_name",
                   attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                   values=SCL_sigName, mart= ensembl_hg38, useCache = FALSE)
SCL_sigID <- subset(SCL_sigID, !startsWith(SCL_sigID$chromosome_name, 'CH'))
SCL_sigID <- subset(SCL_sigID, SCL_sigID$gene_biotype == 'protein_coding')
dim(SCL_sigID)
row.names(SCL_sigID) <- SCL_sigID$external_gene_name

YAP_TAZ_list <- c('CYR61','CTGF','AMOTL2','ANKRD1','IGFBP3','F3','FJX1','NUAK2','LAST2','CRIM1','GADD45A','TGFB2','PTPN14','NTSE','FOXF2','AXL','DOCK5','ASAP1','RBMS3','MYOF','ANHGEF17','CCDC80','AJUBA')
YAP_TAZ_sigID <- getBM(filters= "external_gene_name",
                       attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                       values=YAP_TAZ_list, mart= ensembl_hg38, useCache = FALSE)
YAP_TAZ_sigID <- subset(YAP_TAZ_sigID, !startsWith(YAP_TAZ_sigID$chromosome_name, 'CH'))
row.names(YAP_TAZ_sigID) <- YAP_TAZ_sigID$external_gene_name

basal_list <- c('COL17A1','CSMD2','CDH13','MUM1L1','MMP3','IL33','GIMAP8','PDPN','VSNL1','BNC1','IGFBP7','DLK2','HMGA2','NOTCH4','THBS2','TAGLN','FHL1','ANXA8L2','COL4A6','KCNQ5','WNT7A','KCNMA1','NIPAL4','FLRT2','LTBP2','FOXI1','NGFR','SERPINB13','CNTNAP3B','FGFR3','ARHGAP25','AEBP1','FJX1','TNC','MSRB3','NRG1','SERPINF1','DLC1','IL1A','DKK3','ERG','SYNE1','JAG2','JAM3','MRC2','SPARC','C16orf74','FAT3','KIRREL','SH2D5','KRT6A','KRT34','ITGA6','TP63','KRT5','KRT14')
basal_sigID <- getBM(filters= "external_gene_name",
                     attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                     values=basal_list, mart= ensembl_hg38, useCache = FALSE)
basal_sigID <- subset(basal_sigID, !startsWith(basal_sigID$chromosome_name, 'CH'))
row.names(basal_sigID) <- basal_sigID$external_gene_name

row.names(rlog)[which(row.names(rlog) %in% AR_sigID$external_gene_name)] <- AR_sigID[row.names(rlog)[which(row.names(rlog) %in% AR_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% NEPC_sigID$external_gene_name)] <- NEPC_sigID[row.names(rlog)[which(row.names(rlog) %in% NEPC_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% WNT_sigID$external_gene_name)] <- WNT_sigID[row.names(rlog)[which(row.names(rlog) %in% WNT_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% SCL_sigID$external_gene_name)] <- SCL_sigID[row.names(rlog)[which(row.names(rlog) %in% SCL_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% YAP_TAZ_sigID$external_gene_name)] <- YAP_TAZ_sigID[row.names(rlog)[which(row.names(rlog) %in% YAP_TAZ_sigID$external_gene_name)],'ensembl_gene_id_version']
row.names(rlog)[which(row.names(rlog) %in% basal_sigID$external_gene_name)] <- basal_sigID[row.names(rlog)[which(row.names(rlog) %in% basal_sigID$external_gene_name)],'ensembl_gene_id_version']


gsvaScore <- gsva(as.matrix(rlog), list(AR_sigID$ensembl_gene_id_version, WNT_sigID$ensembl_gene_id_version, NEPC_sigID$ensembl_gene_id_version, SCL_sigID$ensembl_gene_id_version), mx.diff=F, min.sz=5, max.sz=500, method="gsva")
gsvaScore
patient_combat_data_result[,c("AR_score","WNT_score","NEPC_score","SCL_score")] <- t(gsvaScore)[rownames(patient_combat_data_result),]

patient_combat_data_result[,c("AR_score","WNT_score","NEPC_score","SCL_score",'YAP_TAZ_score','basal_score')] <- cbind(colSums(rlog[AR_sigID$ensembl_gene_id_version,]),colSums(rlog[WNT_sigID$ensembl_gene_id_version,]),colSums(rlog[NEPC_sigID$ensembl_gene_id_version,]),colSums(rlog[SCL_sigID$ensembl_gene_id_version,]),colSums(rlog[YAP_TAZ_sigID$ensembl_gene_id_version,]),colSums(rlog[basal_sigID$ensembl_gene_id_version,]))
# write.csv(patient_combat_data_result, '../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_su2cPatients_noPDX_noBatch_wSignalScore.csv')
write.csv(patient_combat_data_result, '../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_su2cPatients_top93_wSignalScore.csv')

##########################################################################################################################################################
##########################################################################################################################################################
# Figure S12a

# Alternative function for simple annotation
col = c("DeepDel" = "blue", "AMP" = "red", "mutation" = "#008000")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  DeepDel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["DeepDel"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["mutation"], col = NA))
  }
)


mydata<-read.csv('../EX_0208_prostateRnaCluster/PNAS_su2c/transformed_su2c.csv',header=TRUE, row.names = 1, check.names=T)
oncoPrint(mydata[,patient_order],alter_fun=alter_fun,col=col,show_column_names = TRUE,
          remove_empty_rows = F, row_order = row.names(mydata), column_order = colnames(mydata))

# patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/ericaSigGene_ntpPrediction_su2cPatients_noPDX_noBatch.csv', row.names = 1, stringsAsFactors = F, check.names=T)
patient_combat_data_result <- read.csv('../EX_0208_prostateRnaCluster/template_byTop/RES_rankBy_baseMean_top93_su2c.csv', row.names = 1, stringsAsFactors = F, check.names=T)
dim(patient_combat_data_result)

patient_combat_data_result[which(patient_combat_data_result$FDR >= 0.5),'prediction'] <- 'unknown'
patient_combat_data_result[which(patient_combat_data_result$p.value >= 0.05),'prediction'] <- 'unknown'


ar_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'AR_dependent'),])
wnt_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'WNT'),])
nepc_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'NEPC'),])
scl_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'SCL'),])
unknown_patient <- row.names(patient_combat_data_result[which(patient_combat_data_result$prediction == 'unknown'),])

patient_order <- c(ar_patient, wnt_patient, nepc_patient, scl_patient,unknown_patient)

dim(mydata[,patient_order])
mydata <- mydata[,patient_order]
oncoPrint(mydata,alter_fun=alter_fun,col=col,show_column_names = TRUE,
          remove_empty_rows = F, row_order = row.names(mydata), column_order = colnames(mydata))

metaData <- read.table('../EX_0208_prostateRnaCluster/PNAS_su2c/data_clinical_sample.csv', sep=',', row.names = 1, header = T, check.names=TRUE)
rownames(metaData) <- make.names(rownames(metaData))
metaData[patient_order,'ETS_FUSION_SEQ']
arScores <- metaData[patient_order,'AR_SCORE']
arScores <- arScores[!is.na(arScores)]

metaData['group'] <- patient_combat_data_result[rownames(metaData),'prediction']

col_AR <- colorRamp2(c(quantile(arScores, 0), quantile(arScores, 1)), c("white", '#F217F7'))
ha = HeatmapAnnotation(df = metaData[patient_order, c('group','AR_SCORE','ETS_FUSION_SEQ')], 
                       col = list(group = c("AR_dependent" = "#F8766D", "WNT" = "#7CAE00", "NEPC" = "#00BFC4", "SCL"= "#C77CFF", 'unknown' = 'grey'), AR_SCORE=col_AR, ETS_FUSION_SEQ=c("Negative"="grey", "Positive"="orange")))

column_title = "OncoPrint"
heatmap_legend_param = list(title = "Alternations", at = c("DeepDel", "AMP", "mutation"), 
                            labels = c("Deep deletion", "Amplification", "Mutation"))
oncoPrint(mydata,
          alter_fun = alter_fun, col = col, 
          top_annotation = ha,
          row_order = row.names(mydata), 
          column_order = patient_order,
          column_title = column_title, 
          heatmap_legend_param = heatmap_legend_param,
          show_column_names = T,
          #remove_empty_rows = T,
          show_heatmap_legend=T,
          left_annotation =  rowAnnotation(
            rbar = anno_oncoprint_barplot(
              axis_param = list(direction = "reverse")
            )),
          right_annotation = NULL)


##########################################################################################################################################################
# Figure S12b
mydata<-read.csv('../EX_0208_prostateRnaCluster/PNAS_su2c/transformed_su2c_addZNRF3.csv',header=TRUE, row.names = 1, check.names=T)

# wnt_data <- mydata[c('APC','CTNNB1','RNF43','ZNRF3','TCF7L2','RSPO2','HDAC2'), wnt_patient]
wnt_data <- mydata[c('APC','CTNNB1','ZNRF3','TCF7L2','RSPO2'), wnt_patient]

oncoPrint(wnt_data,alter_fun=alter_fun,col=col,show_column_names = TRUE,
          remove_empty_rows = F, row_order = row.names(wnt_data), column_order = colnames(wnt_data))


##########################################################################################################################################################
# Figure 4D

fy_plot <- readRDS('../EX_0208_prostateRnaCluster/PNAS_su2c/WNT_MUT_ggplot2.RDS')
fy_plot

metaData <- read.table('../EX_0208_prostateRnaCluster/PNAS_su2c/data_clinical_sample.csv', sep=',', row.names = 1, header = T, check.names=TRUE)
su2c_rNames <- make.names(rownames(metaData))

write.csv(cbind(su2c_rNames,rownames(metaData)),'../EX_0208_prostateRnaCluster/su2c_rName.csv', quote=F)


#######################################################################################################################################################################################################################################################################################
# Figure 4E

sur_su2c <- read.csv('../EX_0208_prostateRnaCluster/PNAS_su2c/su2c_survival_r_top93.csv', row.names = 1, check.names = F)
sur_su2c <- read.csv('../EX_0208_prostateRnaCluster/PNAS_su2c/Table_S11_Patient_survival.csv', row.names = 1, check.names = F)

surv_object_1 <- Surv(time = sur_su2c$OS_MONTHS, event = sur_su2c$OS_STATUS_i)

fit1 <- survfit(surv_object_1 ~ ntpPrediction , data = sur_su2c)
summary(fit1)
ggsurvplot(fit1, 
           data = sur_su2c,
           pval = T)


table(sur_su2c[which(sur_su2c$OS_MONTHS != 'NA'),'ntpPrediction'])
# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
#      30        7       24        3       17

# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
# 39        5       17        3       17

sur_su2c_ar <- sur_su2c[which(sur_su2c$ntpPrediction%in%c("CRPC-AR", "CRPC-SCL")),]

table(sur_su2c_ar[which(sur_su2c_ar$OS_MONTHS != 'NA'),'ntpPrediction'])
# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
# 30        0       24        0        0

# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
# 39        0       17        0        0 

surv_object_2 <- Surv(time = sur_su2c_ar$OS_MONTHS, event = sur_su2c_ar$`OS_STATUS_i`)

fit2 <- survfit(surv_object_2 ~ ntpPrediction , data = sur_su2c_ar)
summary(fit2)
ggsurvplot(fit2, 
           data = sur_su2c_ar,
           pval = T)


table(sur_su2c_ar[which(sur_su2c_ar$OFF_ARSI != 'NA'),'ntpPrediction'])
# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
# 22        0       16        0        0 

# CRPC-AR  CRPC-NE CRPC-SCL CRPC-WNT      UNK 
# 29        0        9        0        0

surv_object_3 <- Surv(time = sur_su2c_ar$`Time on first-line ARSI`, event = sur_su2c_ar$OFF_ARSI)
surv_object_3

fit3 <- survfit(surv_object_3 ~ ntpPrediction , data = sur_su2c_ar)
summary(fit3)

ggsurvplot(fit3, 
           data = sur_su2c_ar,
           pval = T)

dev.cur()
dev.off(1)

##########################################################################################################################################################
# Make a version of su2c with hg19 gene names for SVC model prediction
su2c <- read.csv('../EX_0208_prostateRnaCluster/su2c_rnaseq_geneNameRevised_noPDX_noBatch.csv', header = T, row.names = 1)



# # Use scaled, not centered
# # below is for centered
# ht_list_2 <- Heatmap(as.matrix(mat_scaled_2), 
#                      top_annotation = ha,
#                      cluster_rows = F,
#                      cluster_columns = F,
#                      show_row_names = F,
#                      show_column_names = F,
#                      show_heatmap_legend = F,
#                      column_title = "Molecular classification of IPM patients", 
#                      column_title_gp = gpar(fontsize = 14, fontface = "bold")) + ha_row
# saveRDS(ht_list, paste0("/Users/fanyingtang/Downloads/figures/", "IPM_NTP_results_centered-", Sys.Date(), ".RDS"))
# pdf(paste0("/Users/fanyingtang/Downloads/figures/", "IPM_NTP_results_centered-", Sys.Date(), ".pdf"), width=6, height=6, useDingbats=F)
# draw(ht_list_2, heatmap_legend_side = "right", annotation_legend_side = "left")
# dev.off()


##########################################################################################################################################################
# Archived
patient_combat_data <- read.csv('../EX_0201_prostateSubtypePredict/in_house_patients/data_mRNA_seq_fpkm_polya.txt', header = T, sep='\t')
# patient_combat_data <- matrixData[rowSums(patient_combat_data) > 0, ]
geneList <- getBM(filters= "external_gene_name",
                  attributes= c("ensembl_gene_id_version","ensembl_gene_id","chromosome_name", "external_gene_name", "gene_biotype"),
                  values=patient_combat_data$Hugo_symbol, mart= ensembl, useCache = FALSE)
geneList <- subset(geneList, !startsWith(geneList$chromosome_name, 'H'))
geneList <- subset(geneList, !startsWith(geneList$chromosome_name, 'GL'))
geneList <- unique(geneList)
geneList <- subset(geneList, geneList$gene_biotype == "protein_coding")
patient_combat_data <- patient_combat_data[patient_combat_data$Hugo_symbol %in% geneList$external_gene_name,]
patient_combat_data <- unique(patient_combat_data)
row.names(patient_combat_data) <- patient_combat_data$Hugo_symbol

repeat_gene_name <- c("AKAP17A","ASMT","ASMTL","CD99","CRLF2","CSF2RA","DHRSX","GTPBP6","IL3RA","IL9R","P2RY8","PLCXD1","PPP2R3B","SHOX","SLC25A6","SPRY3","VAMP7","ZBED1")
removeList <- c()
for(x in which(patient_combat_data$Hugo_symbol %in% repeat_gene_name)){
  if(sum(as.numeric(patient_combat_data[x,2:271])) == 0){
    removeList <- c(removeList,x)
  }
}
removeList
patient_combat_data <- patient_combat_data[-removeList, ]
row.names(patient_combat_data) <- patient_combat_data$Hugo_symbol

