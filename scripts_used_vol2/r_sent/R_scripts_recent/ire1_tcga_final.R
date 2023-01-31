#Rscript /Users/user/tcga_ire1.R /Users/user/GDCdata/RData/TCGA-BRCA.RData /Users/user/TCGA-BRCA.txt
rm(list = ls()) 

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(getopt))

args = commandArgs(trailingOnly=TRUE)



load(args[1])
#load("/Users/user/GDCdata/RData/TCGA-BRCA.RData")
# a<-"/Users/user/GDCdata/RData/TCGA-BRCA.RData"
# 
# a_split<-strsplit(a,"/")[[1]]
# b<-a_split[length(a_split)]
# b<-strsplit(b,"\.")

# xbp1<-read.table("/Users/user/Downloads/signatures_top_secret/xbp1.csv", sep=",")
# FC<-ifelse(xbp1$V3=="TRUE", as.numeric(1), as.numeric(-1))
# xbp1<-as.data.frame(cbind(xbp1$V2,FC))
# 
# ridd<-read.table("/Users/user/Downloads/signatures_top_secret/ridd.csv", sep=",")
# FC<-ifelse(ridd$V3=="TRUE", as.numeric(1), as.numeric(-1))
# ridd<-as.data.frame(cbind(ridd$V2,FC))
# 
# genesignx<-xbp1[,1]
# genesignr<-ridd[,1]
# genesignx.FC<-xbp1
# genesignr.FC<-ridd


sign_xbp1=c("ASS1", "C3", "CCL20", "COL4A6", "CXCL2", "CXCL5", "CXCL8", "IFI44L", "IL1B", "IL6", "KCNN2", "MMP1", "MMP12", "MMP3", "PLA2G4A", "PPP4R4", "SERPINB2", "TFPI2",
            "ZNF804A")
genesignx.FC<-as.data.frame(cbind(sign_xbp1, rep(as.numeric(1), length(sign_xbp1))))
genesignx<-sign_xbp1
sign_ridd=c("ANGPT1", "CFH", "CFI", "CLEC3B", "COL3A1", "COL8A1", "DACH1", "DCN", "FHL1", "GAS1", "LUM", "OXTR", "PLAC8", "RGS4", "TAGLN", "TGFB2", "THBS1", "TIMP3", "TMEM255A")
genesignr<-sign_ridd
genesignr.FC<-as.data.frame(cbind(sign_ridd, rep(as.numeric(-1), length(sign_xbp1))))



not.found <- genesignx[!(genesignx%in% rownames(exp.filt))]
if(length(not.found) > 0){
  msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
  print(msg)
}

not.found <- genesignr[!(genesignr%in% rownames(exp.filt))]
if(length(not.found) > 0){
  msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
  print(msg)
}

#keep the genes of signature that are present in exp.filt for further analysis
genesignx <-  genesignx[genesignx %in% rownames(logCPM)]
genesignr<-  genesignr[genesignr %in% rownames(logCPM)]


#map the genesign on logCPM.matrix
logCPM.genesignx <- logCPM[genesignx, ]
logCPM.genesignr <- logCPM[genesignr, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
qx <- rowQuantiles(logCPM.genesignx, probs=probs)
qr <- rowQuantiles(logCPM.genesignr, probs=probs)
#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
auxx <- matrix(rep(0), nrow=nrow(logCPM.genesignx), ncol=ncol(logCPM.genesignx))
rownames(auxx) <- rownames(logCPM.genesignx)
colnames(auxx) <- colnames(logCPM.genesignx)
#fill in the cells of aux matrix
for(j in 1:ncol(auxx)){
  for(i in 1:nrow(auxx)){
    if(genesignx.FC[i,2] > 0) {   
      if (logCPM.genesignx[i,j] <= qx[i,1]){ 
        auxx[i,j] <- 1
      }else if (logCPM.genesignx[i,j] > qx[i,1] & logCPM.genesignx[i,j] <= qx[i,2]){
        auxx[i,j] <- 2
      }else if(logCPM.genesignx[i,j] > qx[i,2] & logCPM.genesignx[i,j] < qx[i,3]){
        auxx[i,j] <- 3
      }else auxx[i,j] <- 4
    }
    else if(genesignx.FC[i,2] < 0){
      if (logCPM.genesignx[i,j] <= qx[i,1]){ 
        auxx[i,j] <- 4
      }else if (logCPM.genesignx[i,j] > qx[i,1] & logCPM.genesignx[i,j] <= qx[i,2]){
        auxx[i,j] <- 3
      }else if(logCPM.genesignx[i,j] > qx[i,2] & logCPM.genesignx[i,j] < qx[i,3]){
        auxx[i,j] <- 2
      }else auxx[i,j] <- 1
    }
  }
}

scorex <- colSums(auxx)/length(genesignx)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

auxr <- matrix(rep(0), nrow=nrow(logCPM.genesignr), ncol=ncol(logCPM.genesignr))
rownames(auxr) <- rownames(logCPM.genesignr)
colnames(auxr) <- colnames(logCPM.genesignr)
#fill in the cells of aux matrix
for(j in 1:ncol(auxr)){
  for(i in 1:nrow(auxr)){
    if(genesignr.FC[i,2] > 0) {   
      if (logCPM.genesignr[i,j] <= qr[i,1]){ 
        auxr[i,j] <- 1
      }else if (logCPM.genesignr[i,j] > qr[i,1] & logCPM.genesignr[i,j] <= qr[i,2]){
        auxr[i,j] <- 2
      }else if(logCPM.genesignr[i,j] > qr[i,2] & logCPM.genesignr[i,j] < qr[i,3]){
        auxr[i,j] <- 3
      }else auxr[i,j] <- 4
    }
    else if(genesignr.FC[i,2] < 0){
      if (logCPM.genesignr[i,j] <= qr[i,1]){ 
        auxr[i,j] <- 4
      }else if (logCPM.genesignr[i,j] > qr[i,1] & logCPM.genesignr[i,j] <= qr[i,2]){
        auxr[i,j] <- 3
      }else if(logCPM.genesignr[i,j] > qr[i,2] & logCPM.genesignr[i,j] < qr[i,3]){
        auxr[i,j] <- 2
      }else auxr[i,j] <- 1
    }
  }
}
#calculate a patient score
scorer <- colSums(auxr)/length(genesignr)

colData(exp)$level.group.xbp1 <- "medium"
min.cut <- max(sort(scorex)[1:(length(scorex) * threshold)])
high.cut <- min(sort(scorex, decreasing = T)[1:(length(scorex) * threshold)])
colData(exp)[scorex <= min.cut,"level.group.xbp1"] <- "xbp1-"
colData(exp)[scorex >= high.cut,"level.group.xbp1"] <- "xbp1+"


colData(exp)$level.group.ridd <- "medium"
min.cut <- max(sort(scorer)[1:(length(scorer) * threshold)])
high.cut <- min(sort(scorer, decreasing = T)[1:(length(scorer) * threshold)])
colData(exp)[scorer <= min.cut,"level.group.ridd"] <- "ridd+"
colData(exp)[scorer >= high.cut,"level.group.ridd"] <- "ridd-"


colData(exp)$grouping<-paste(colData(exp)$level.group.xbp1, colData(exp)$level.group.ridd, sep="_")  #!!!!!!!
colData(exp)$grouping2<-colData(exp)$grouping

test1<-as.data.frame(cbind(colData(exp)$barcode,colData(exp)$grouping2))
index1<-which(test1$V2=="xbp1+_ridd+")
index2<-which(test1$V2=="xbp1+_ridd-")
index3<-which(test1$V2=="xbp1-_ridd-")
index4<-which(test1$V2=="xbp1-_ridd+")
test2<-rbind(test1[index1,],test1[index2,],test1[index3,],test1[index4,])
test2$V3<-c(rep("A", length(index1)), rep("B", length(index2)), rep("C", length(index3)), rep("D", length(index4)))
write.table(test2, args[2], sep="\t", quote=FALSE)
