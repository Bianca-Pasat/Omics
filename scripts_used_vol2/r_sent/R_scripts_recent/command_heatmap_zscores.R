
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler)) #
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(circlize)) 
suppressPackageStartupMessages(library(ComplexHeatmap)) #
suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(getopt))

logCPM.genesign <- expHif

zscores1<-scale(exp)
test<-zscores1[rownames(zscores1) %in% rownames(expHif),]
zscores2<-scale(logCPM.genesign)
test<-zscores2

col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5

clinical_patient_Cancer<-clinical_patient_Cancer[rownames(clinical_patient_Cancer) %in% colnames(test),]

signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(sample_type=sample_info$condition)

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
pdf("heatOutPdf")
hm<-Heatmap(test, name="expression",
        col=col_fun,
        cluster_rows=T,
        cluster_columns = F,
        row_names_side = "left",
        show_row_names = T,
        #column_order = order(clinical_patient_Cancer$shortLetterCode),
        #row_order = order(rownames(logCPM.genesign)),
        show_column_names = F,
        #column_names_gp=gpar(cex=fontsize),
        show_row_dend = FALSE,
        show_column_dend=TRUE,
        row_names_gp=gpar(cex=fontsize+0.1),
        row_names_max_width = unit(5, "cm"),
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "ward.D",
        clustering_distance_columns =  "euclidean",
        clustering_method_columns = "ward.D",
        row_dend_width = unit(10, "mm"),
        #row_km = 2,
        #column_km =2,
        #km=2,
        left_annotation=row,
        bottom_annotation = col,
        heatmap_width = unit(10, "cm"),
        width = NULL,
        heatmap_height = unit(10, "cm"))
draw(hm)
