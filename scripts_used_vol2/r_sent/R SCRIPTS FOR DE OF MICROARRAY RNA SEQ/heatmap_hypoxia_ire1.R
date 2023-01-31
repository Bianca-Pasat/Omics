# HYPOXIA AND MKC
suppressPackageStartupMessages(library(clusterProfiler)) #
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(circlize)) 
suppressPackageStartupMessages(library(ComplexHeatmap)) #
suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(getopt))


hif<-read.table("//home/bianca/Desktop/hypoxia-claire/hif_signature.txt",sep="\t",row.names=1)
hypoxia<-read.table("//home/bianca/Desktop/hypoxia-claire/dmso_hyp_4_prioritized.tsv",sep="\t",row.names=1)
hypoxia_all<-read.table("//home/bianca/Desktop/hypoxia-claire/hypoxia_dmso_4_DEregulatedBio_NEW.txt",sep="\t",row.names=1)
ire1_all<-read.table("/home/bianca/Downloads/ANALYSIS_APRIL/DE/MICROARRAY/microarray_rankprod_DE_Bio_4.txt",sep="\t",row.names=1)
ire1_sig<-read.table("/home/bianca/Downloads/ANALYSIS_APRIL/DE/MICROARRAY/dmso_vs_mkc_4_de.tsv",sep="\t",row.names=1)

library(biomaRt)
httr::set_config(httr::config(ssl_verifypeer = FALSE))



homo.anno = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")


values<-rownames(ire1_all)
attributes <- c("ensembl_gene_id", "hgnc_symbol")


filters="ensembl_gene_id"

# #listAttributes(homo.anno)
genes_ire1 <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
ire1_all$genes=rownames(ire1_all)
genes_full_ire1=merge(ire1_all, genes_ire1, by.x=3, by.y=1 )

hypoxia$genes=rownames(hypoxia)
hypoxia_and_ire1<-merge(genes_full_ire1, hypoxia, by.x=4,by.y=1)
hypoxia_and_ire1<-as.data.frame(hypoxia_and_ire1[,c(1,3,4,9,10)])
colnames(hypoxia_and_ire1)<-c("Gene_name","FC_MKC","pval_MKC","FC_Hypoxia","pval_Hypoxia")
rownames(hypoxia_and_ire1)<-hypoxia_and_ire1$Gene_name

hypoxia_all$genes=rownames(hypoxia_all)
hypoxia_and_ire1<-merge(genes_full_ire1, hypoxia_all, by.x=4,by.y=3)
hypoxia_and_ire1<-hypoxia_and_ire1[,-2]
colnames(hypoxia_and_ire1)<-c("Gene_name","FC_MKC","pval_MKC","FC_Hypoxia","pval_Hypoxia")
rownames(hypoxia_and_ire1)<-hypoxia_and_ire1$Gene_name


hypoxia_all$genes<-rownames(hypoxia_all)
hyp_ire<-merge(hypoxia_all, ire1_sig, by.x=3,by.y=1)
hyp_ire=hyp_ire[,c(1,2,8)]
hyp_ire$V7<-as.numeric(hyp_ire$V7)
rownames(hyp_ire)<-hyp_ire$genes
hyp_ire=hyp_ire[,-1]
colnames(hyp_ire)<-c("FC_Hypoxia","FC_MKC")
test<-hyp_ire
test$FC_Hypoxia_x10=test$FC_Hypoxia * 10
test=test[,c(2,3)]

## significant ire1 and FC
test=hypoxia_and_ire1[rownames(hypoxia_and_ire1) %in% ire1_sig$V2,]
test=as.matrix(test[,c(2,4)])
test$FC_Hypoxia_x10=test$FC_Hypoxia*10
test=as.matrix(test[,c(2,6)])

## significant hypoxia and FC
test=hypoxia_and_ire1[rownames(hypoxia_and_ire1) %in% hypoxia$V2,]
test=as.matrix(test[,c(2,4)])
test$FC_Hypoxia_x10=test$FC_Hypoxia*10
test=as.matrix(test[,c(2,6)])

## all common genes
test=hypoxia_and_ire1
test=as.matrix(test[,c(2,4)])
test$FC_Hypoxia_x10=test$FC_Hypoxia*10
test=as.matrix(test[,c(2,6)])

### proteomics
test=as.matrix(de_all)
d8=read.table("/home/bianca/Desktop/deModules_with_propr_and_rankprod_8_FC.txt", sep="\t")
d24=read.table("/home/bianca/Desktop/deModules_with_propr_and_rankprod_24_FC.txt", sep="\t")
d48=read.table("/home/bianca/Desktop/deModules_with_propr_and_rankprod_48_FC.txt", sep="\t")
d72=read.table("/home/bianca/Desktop/deModules_with_propr_and_rankprod_72_FC.txt", sep="\t")

d8$name=rownames(d8)
d24$name=rownames(d24)
d48$name=rownames(d48)
d72$name=rownames(d72)
one<-merge(d8,d24 ,by ="name")
two<-merge(one, d48, by="name")
test<-merge(two, d72, by="name")
rownames(test)=test$name
test<-test[,c(2,4,6,8)]
colnames(test)<-c(8,24,48,72)



col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5


clinical_patient_Cancer<-colnames(test)

signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(Treatment=clinical_patient_Cancer)

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
#pdf(heatOutPdf)
#pdf("sig_ire1_and_FC")
Heatmap(test, name="expression",
        col=col_fun,
        cluster_rows=T,
        cluster_columns = F,
        row_names_side = "left",
        show_row_names = T,
        column_order = colnames(test),
        #row_order = order(rownames(logCPM.genesign)),
        show_column_names = T,
        #column_names_gp=gpar(cex=fontsize),
        show_row_dend = TRUE,
        #show_column_dend=TRUE,
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
        #left_annotation=row,
        bottom_annotation = col,
        #heatmap_width = unit(25, "cm"),
        width = NULL)
        #heatmap_height = unit(15, "cm"))
dev.off()
draw(hm)
