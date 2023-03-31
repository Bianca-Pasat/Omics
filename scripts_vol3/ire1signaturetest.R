
library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(ComplexHeatmap)
library(RankProd)

query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")





GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)
#input SE,genesign,typesample and threshold
se <- data

# set your sample type e.g. typesample = "TP"
typesample <- c("TP", "NT") # (solid tumors)

exp <-se[, TCGAquery_SampleTypes(colnames(se), typesample=typesample)]
#Filtering SE for the user-defined sample type

expTP <- TCGAquery_SampleTypes(barcode=colnames(exp), typesample = "TP")

expNT <- TCGAquery_SampleTypes(barcode=colnames(exp),typesample = "NT")

rNames<-as.data.frame(rownames(exp))
rNames$names<-gsub("\\..*","",rownames(exp))
eg = as.data.frame(bitr(rNames$names,
                        fromType="ENSEMBL",
                        toType="SYMBOL",
                        OrgDb="org.Hs.eg.db")) #drop=TRUE



eg <- eg[!duplicated(eg$SYMBOL),]#unique values
eg <- eg[!is.na(eg$SYMBOL),]#remove NA


eg.ensembl <- as.character(eg$ENSEMBL)
dd1 <- match(eg.ensembl,rNames$names)
exp <- exp[dd1, ]
rownames(exp) <- as.character(eg$SYMBOL)
# dd <- match(eg.ensembl,rownames(exp))
# 
# exp <- exp[dd, ]
# 
# rownames(exp) <- as.character(eg$SYMBOL)






# Filtering for low counts
exp.filt <- TCGAanalyze_Filtering(tabDF = assay(exp),
                                  method = "quantile",qnt.cut =  0.25)

#create DGEList object (edgeR)
y.exp <- DGEList(counts = exp.filt)

#TMM normalization
y.exp <- calcNormFactors(y.exp, method = "TMM")
###LOGCPM
logCPM <- cpm(y.exp, prior.count=5, log=TRUE)


sign_xbp1=c("ASS1", "C3", "CCL20", "COL4A6", "CXCL2", "CXCL5", "CXCL8", "IFI44L", "IL1B", "IL6", "KCNN2", "MMP1", "MMP12", "MMP3", "PLA2G4A", "PPP4R4", "SERPINB2", "TFPI2",
            "ZNF804A")
genesign<-sign_xbp1
sign_ridd=c("ANGPT1", "CFH", "CFI", "CLEC3B", "COL3A1", "COL8A1", "DACH1", "DCN", "FHL1", "GAS1", "LUM", "OXTR", "PLAC8", "RGS4", "TAGLN", "TGFB2", "THBS1", "TIMP3", "TMEM255A")
genesign<-sign_ridd
genesign<-c(sign_xbp1, sign_ridd)


not.found <- genesign[!(genesign%in% rownames(exp.filt))]
if(length(not.found) > 0){
  msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
  print(msg)
}

#keep the genes of signature that are present in exp.filt for further analysis
genesign <-  genesign[genesign %in% rownames(logCPM)]
#genesign.FC <- genesign.FC[rownames(genesign.FC) %in% genesign, ]


#indexlogCPM

#map the genesign on logCPM.matrix

tlogCPM.genesignindex <- which(colnames(t(logCPM)) %in% genesign)

#map the genesign on logCPM.matrix
logCPM.genesign <- logCPM[genesign, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
q <- rowQuantiles(logCPM.genesign, probs=probs)

#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
rownames(aux) <- rownames(logCPM.genesign)
colnames(aux) <- colnames(logCPM.genesign)
#fill in the cells of aux matrix
for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(genesign[i] > 0) {   
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 1
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 2
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 3
      }else aux[i,j] <- 4
    }
    else if(genesign[i] < 0){
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 4
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 3
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 2
      }else aux[i,j] <- 1
    }
  }
}

#calculate a patient score
score <- colSums(aux)/length(genesign)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group <- "Mid expression"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group"] <- "Low expression"
colData(exp)[score >= high.cut,"level.group"] <- "High expression"


clusterCol <- "level.group"
surv.data <- colData(exp)[(colData(exp$level.group) %in% "Low expression"),]
#survival Plot
TCGAanalyze_survival(data = colData(exp),#data = surv.data,
                     clusterCol = clusterCol,
                     filename = "/home/ben/Desktop/xbp1_second.png",#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)


col_fun = colorRamp2(seq(min(aux), max(aux), length = 3), c("blue", "white", "red"))
fontsize=0.5
clinical_patient_Cancer<-colData(se)

# signa<-data.frame(rep("xbp1",19))
# signa<-data.frame(rep("ridd",19))
signa<-data.frame(c(rep("xbp1",19), rep("ridd",19)))
col<-HeatmapAnnotation(df=clinical_patient_Cancer[-1,c(4,58,59)])

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
hm<-Heatmap(aux, name="expression",
            col=col_fun,
            cluster_rows=F,
            cluster_columns = T,
            #row_names_side = "left",
            show_row_names = F,
            #column_order = order(clinical_patient_Cancer$vital_status),
            #row_order = order(rownames(logCPM.genesign)),
            show_column_names = F,
            #column_names_gp=gpar(cex=fontsize),
            show_row_dend = FALSE,
            show_column_dend=TRUE,
            row_names_gp=gpar(cex=fontsize+0.1),
            row_names_max_width = unit(5, "cm"),
            clustering_distance_rows ="pearson",
            clustering_method_rows = "ward.D",
            clustering_distance_columns =  "pearson",
            clustering_method_columns = "ward.D",
            row_dend_width = unit(10, "mm"),
            row_km = 4,
            #column_km =4,
            #km=2,
            left_annotation=row,
            bottom_annotation = col,
            heatmap_width = unit(15, "cm"),
            width = NULL,
            heatmap_height = unit(15, "cm"))

taux<-t(aux)
row<-HeatmapAnnotation(signature=signa[,1] )
HM<-Heatmap(taux, name="expression",
        col=col_fun,
        cluster_rows=T,
        cluster_columns = F,
        #row_names_side = "left",
        show_row_names = F,
        #column_order = order(clinical_patient_Cancer$vital_status),
        #row_order = order(rownames(logCPM.genesign)),
        show_column_names = F,
        #column_names_gp=gpar(cex=fontsize),
        show_row_dend = FALSE,
        show_column_dend=TRUE,
        row_names_gp=gpar(cex=fontsize+0.1),
        row_names_max_width = unit(5, "cm"),
        clustering_distance_rows ="pearson",
        clustering_method_rows = "ward.D",
        clustering_distance_columns =  "pearson",
        clustering_method_columns = "ward.D",
        row_dend_width = unit(10, "mm"),
        row_km = 4,
        #column_km =2,
        #km=2,
        #left_annotation=col,
        bottom_annotation = row,
        heatmap_width = unit(15, "cm"),
        width = NULL,
        heatmap_height = unit(15, "cm"))

r.dend <- row_dend(HM)  #Extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)

clusters<-lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
clusters<-lapply(rcl.list, function(x) length(x))  

# loop to extract genes for each cluster.
 for (i in 1:length(row_order(HM))){
 if (i == 1) {
clu <- t(t(row.names(taux[row_order(HM)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
       colnames(out) <- c("GeneID", "Cluster")
       } else {
        clu <- t(t(row.names(taux[row_order(HM)[[i]],])))
         clu <- cbind(clu, paste("cluster", i, sep=""))
        out <- rbind(out, clu)
        }
 }


clusterone<-data.frame(out[c(1:51),1])
rownames(clusterone)<-clusterone[,1]
clustertwo<-data.frame(out[c(52:96),])
rownames(clustertwo)<-clustertwo[,1]
clusterthree<-data.frame(out[c(97:126),])
rownames(clusterthree)<-clusterthree[,1]
clusterfour<-data.frame(out[c(125:170),])
rownames(clusterfour)<-clusterfour[,1]
xhrh<-aux[,(colnames(aux) %in% rownames(clusterone) )]
#colnames(xhrh)<-rep("xhrh", 51)
xhrl<-aux[,(colnames(aux) %in% rownames(clustertwo) )]
#colnames(xhrl)<-rep("xhrl", 45)
xlrh<-aux[,(colnames(aux) %in% rownames(clusterthree) )]
#colnames(xlrh)<-rep("xlrh", 30)
xlrl<-aux[,(colnames(aux) %in% rownames(clusterfour) )]
#colnames(xlrl)<-rep("xlrl", 46)
colData(exp)$exp.data<- ifelse(colnames(exp) %in% colnames(all), "1", "notreported")


colData(exp)[ colnames(exp) %in% colnames(xhrh),"exp.data"] <- "xhrh"

colData(exp)[ colnames(exp) %in% colnames(xhrl),"exp.data"] <- "xhrl"

colData(exp)[ colnames(exp) %in% colnames(xlrh),"exp.data"] <- "xlrh"

colData(exp)[ colnames(exp) %in% colnames(xlrl),"exp.data"] <- "xlrl"


clusterCol<-"exp.data"
# exp2<-exp[colnames(exp) %in% colnames(all),]
# colData(exp)$expgroup<-low.all
# all<-cbind(xhrh,xhrl,xlrh,xlrl )

TCGAanalyze_survival(data = colData(exp),#data = surv.data,
                     clusterCol = clusterCol,
                     filename = "/Users/user/Desktop/alltogetherLdsurvival.png",#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)




write.table(out, file= "gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)
