###############This code is created by Konstantinos Voutetakis###############################################
#################################SurvivalPlot################################################################

# This module takes as input four arguments: a) a summarizedExperiment (SE) with TCGA gene expression and
# sample data with the following columns included:  days_to_death, vital_status, days_to_last_follow_up,   
# b) a gene signature (genesign) with gene_ids and fold changes, c) a survival threshold 
#############################################################################################################
###Rscript /Users/user/survivalTCGA.R "TCGA-COAD" /Users/user/Downloads/HCT_116_MRNA-MIRNA/survival/network_common/wptb_common_network/Gene_Prioritization.tsv /Users/user/Downloads/HCT_116_MRNA-MIRNA/survival/network_common/wptb_common_network.txt /Users/user/wptb_network_survival.jpg

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
args = commandArgs(trailingOnly=TRUE)


########## bopmart
mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "uswest")
attributes = c("ensembl_gene_id", "hgnc_symbol")
filters="hgnc_symbol"
################
XBP1sign <- c("ASS1","C3","CCL20","COL4A6","CXCL2","CXCL5","CXCL8","IFI44L","IL1B","IL6","KCNN2","MMP1","MMP12","MMP3","PLA2G4A","PPP4R4","SERPINB2","TFPI2","ZNF804A")
RIDDsign <- c("ANGPT1","CFH","CFI","CLEC3B","COL3A1","COL8A1","DACH1","DCN","FHL1","GAS1","LUM","OXTR","PLAC8","RGS4","TAGLN","TGFB2","THBS1","TIMP3","TMEM255A")

values=rownames(exp)
expENS <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

# write.table(XBP1sign, "/home/ben/Desktop/xbp1.txt", sep="\t", quote = F, row.names = F)
# write.table(RIDDsign, "/home/ben/Desktop/ridd.txt", sep="\t", quote = F, row.names = F)
XBP1sign<-read.table("/home/ben/Desktop/xbp1.txt", sep = "\t")
RIDDsign<-read.table("/home/ben/Desktop/ridd.txt", sep = "\t")

#input SE,genesign,typesample and threshold
p<-paste(text=args[1])

## see wehre you can keep only basal
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")


GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)
#input SE,genesign,typesample and threshold
se <- data

# set your sample type e.g. typesample = "TP"
typesample <- c("TP", "NT") # (solid tumors)
typesample <- c("TP") # (solid tumors)
all<-colData(se)
tnbc<-all[which(all$paper_BRCA_Subtype_PAM50=="Basal"),]
se1<-se[,colnames(se) %in% tnbc$barcode]
exp <-se1[, TCGAquery_SampleTypes(colnames(se1), typesample=typesample)]

#Filtering SE for the user-defined sample type

expTP <- TCGAquery_SampleTypes(barcode=colnames(exp), typesample = "TP")

expNT <- TCGAquery_SampleTypes(barcode=colnames(exp),typesample = "NT")

genesign<-read.table(args[2],  header = T,  sep="\t")[,2]
genesign.FC<-read.table(args[3],  header = F, row.names=1, sep="\t")[,c(1,2)]

# set your sample type e.g. typesample = "TP"
#typesample <- "TP" # (solid tumors)


#Filtering SE for the user-defined sample type
#exp <- se[, TCGAquery_SampleTypes(colnames(se), typesample=typesample)]

# #map and annotate gene symbols
# #ensembl gene ids to external gene symbols
# 
# eg = as.data.frame(bitr(rownames(genesign.FC),
#                         fromType="ENSEMBL",
#                         toType="SYMBOL",
#                         OrgDb="org.Hs.eg.db"))
# 
# eg <- eg[!duplicated(eg$SYMBOL),]#unique values
# eg <- eg[!is.na(eg$SYMBOL),]#remove NA
# 
# 
# eg.ensembl <- as.character(eg$ENSEMBL)
# dd <- match(eg.ensembl,rownames(genesign.FC))
# 
# genesign.FC <-genesign.FC[dd ,]
# rownames(genesign.FC) <- as.character(eg$SYMBOL)
# genesign.FC2<-as.data.frame(genesign.FC[,-2])
# rownames(genesign.FC2)<-rownames(genesign.FC)
# genesign.FC<-genesign.FC2
# 
# eg = as.data.frame(bitr(rownames(exp),
#                         fromType="ENSEMBL",
#                         toType="SYMBOL",
#                         OrgDb="org.Hs.eg.db")) #drop=TRUE
# 
# eg <- eg[!duplicated(eg$SYMBOL),]#unique values
# eg <- eg[!is.na(eg$SYMBOL),]#remove NA 
# 
# 
# eg.ensembl <- as.character(eg$ENSEMBL)
# dd <- match(eg.ensembl,rownames(exp))
# 
# exp <- exp[dd, ]
# rownames(exp) <- as.character(eg$SYMBOL)

# Filtering for low counts
exp.filt <- TCGAanalyze_Filtering(tabDF = assay(exp),
                                  method = "quantile",qnt.cut =  0.3)


#Print a message for genes of genesign were not found in exp.filt
not.found <- genesign[!(genesign %in% rownames(exp.filt))]
if(length(not.found) > 0){
  msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
  print(msg)
}

#keep the genes of signature that are present in exp.filt for further analysis
genesign <-  genesign[genesign %in% rownames(exp.filt)]
genesign.FC <- genesign.FC[rownames(genesign.FC) %in% genesign,]

#create DGEList object (edgeR)
y.exp <- DGEList(counts = exp.filt)

#TMM normalization
y.exp <- calcNormFactors(y.exp, method = "TMM")

#calculate log2CMP values
logCPM <- cpm(y.exp, prior.count=5, log=TRUE)

# remove .number from rownames
vv<-gsub("\\..*", "", rownames(logCPM))
rownames(logCPM)<-vv


#map the genesign on logCPM.matrix
genesign<-c(XBP1sign$V1, RIDDsign$V1)
logCPM.genesign <- as.data.frame(logCPM[rownames(logCPM) %in% genesign, ])


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
q <- rowQuantiles(logCPM, probs=probs)
#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
# aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
# rownames(aux) <- rownames(logCPM.genesign)
# colnames(aux) <- colnames(logCPM.genesign)
aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
rownames(aux) <- rownames(logCPM.genesign)
colnames(aux) <- colnames(logCPM.genesign)
#fill in the cells of aux matrix
for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(rownames(aux)[i] %in% XBP1sign$V1){
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 1
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 2
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 3
      }else if (logCPM.genesign[i,j] >= q[i,3]){
        aux[i,j] <- 4
      }
    }else if(rownames(aux)[i] %in% RIDDsign$V1){
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 4
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 3
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 2
      }else if (logCPM.genesign[i,j] >= q[i,3]){
        aux[i,j] <- 1
      }
    }else {
      next
    }
  }
}

#calculate a patient score
score <- colSums(aux)/(length(genesign) -2)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group <- "Mid expression"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group"] <- "Low expression"
colData(exp)[score >= high.cut,"level.group"] <- "High expression"




surv.data <- colData(exp)[colData(exp)$level.group %in% c("Low expression", "High expression"),]
clusterCol <- "level.group"




# d<-as.data.frame(cbind(colData(exp)$vital_status,colData(exp)$days_to_death,colData(exp)$days_to_last_follow_up))
# rownames(d)<-colData(exp)$barcode
# d<-na.omit(d)
# colnames(d)<-c('vital_status','days_to_death','days_to_last_follow_up')
# complete_death<-match(rownames(d),colData(exp)$barcode)


####THIS IS HOW YOU CHOOSE WHAT YOU WANT FROM EXP
# surv.data <- colData(exp)[colData(exp)$barcode %in% rownames(d),]

#exp <-exp[complete_death]
#survival Plot
TCGAanalyze_survival(data = colData(exp),
                     clusterCol = clusterCol,
                     filename = args[4],#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)

TCGAanalyze_survival(data = colData(exp),
                     clusterCol = clusterCol,
                    #filename = args[4],#or "survivalPlot.png"s
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = FALSE,
                     risk.table = TRUE,
                     conf.int= FALSE)


####HEATMAP

colData(exp)$heat <- as.numeric("0")
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"heat"] <- as.numeric("-1")
colData(exp)[score >= high.cut,"heat"] <- as.numeric("1")


#zscores<-scale(logcpmgenexr)


testx<-matrix(colData(exp)$heat, nrow=nrow(exp), ncol=length(colData(exp)$barcode))
rownames(testx) <- rownames(exp)
colnames(testx) <- colnames(logCPM.genesignx)
testx<-testx[genesignx,]

test<-testx


col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5
clinical_patient_Cancer<-colData(se)
clinical_patient_Cancer<-clinical_patient_Cancer[rownames(clinical_patient_Cancer) %in% colnames(test),]

signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(df=clinical_patient_Cancer[,c(4,52)])

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
hm<-Heatmap(test, name="expression",
            col=col_fun,
            cluster_rows=T,
            cluster_columns = T,
            #row_names_side = "left",
            show_row_names = F,
            column_order = order(clinical_patient_Cancer$shortLetterCode),
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
            row_km = 2,
            column_km =2,
            #km=2,
            left_annotation=row,
            bottom_annotation = col,
            heatmap_width = unit(20, "cm"),
            width = NULL,
            heatmap_height = unit(20, "cm"))
draw(hm)
