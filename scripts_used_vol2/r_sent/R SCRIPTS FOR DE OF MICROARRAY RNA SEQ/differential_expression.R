library(tidyverse)
library(GEOquery)
library(preprocessCore)
library(RankProd)
library(ggfortify)
library(ggplot2)
library(sva)
library(biomaRt)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(hugene20sttranscriptcluster.db) # logically the annotation database to query
library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library("RColorBrewer")
library("ggplot2")
library("stringr")
library("reshape2")
library("pheatmap")

##### DIFFERENTIAL EXPRESSION WITH RANKPROD
## rnaseq 8 hours

samples=adjusted
cl<-samples_info$condition

######### 8 HOURS BATCH CORRECTION ON ALL OF THEM
samples=adjusted[,c(1,4,5,8,9,10)]
cl<-samples_info$condition[c(1,4,5,8,9,10)]

######### 24 HOURS BATCH CORRECTION ON ALL OF THEM
samples=adjusted[,c(2,3,6,7,11,12)]
cl<-samples_info$condition[c(2,3,6,7,11,12)]


### miRNA 8 hours
samples=adjusted[,c(-3,-6)]
cl<-samples_info$condition[c(-3,-6)]

#microarray hypoxia
#normalized and removed low count
samples=y$counts
cl<-samples_info$lim2
# hypoxia at dmso
samples<-y$counts[,c(1,2,3,7,8,9)]
cl<-samples_info$lim1[c(1,2,3,7,8,9)]
# hypoxia and mkc
samples<-y$counts[,c(4,5,6,10,11,12)]
cl<-samples_info$lim[c(4,5,6,10,11,12)]

#microarray
samples=adjusted
cl<-samples_vector_file$condition


### 24 hours
samples=adjusted
cl<-samples_vector_file$condition[c(1:3,7:9)]

# rna seq sample info
sample_info_RNAseq$rankprod<-c(rep(1,6), rep(0,6))

#rna seq mrna 8 h
samples=log_cleaned_count[,c(1,4,5,8,9,10)]
cl<-sample_info_RNAseq$rankprod[c(1,4,5,8,9,10)]

#rna seq mrna 24 h
samples=log_cleaned_count[,c(2,3,6,7,11,12)]
cl<-sample_info_RNAseq$rankprod[c(2,3,6,7,11,12)]

## proteomics 8
samples=norm_proteins[,c(10:12,22:24)]
cl=c(rep(1,3),rep(0,3))
# proteomics 8 hours
out_dir="/home/bianca/Desktop/rp_prot_mkc_vs_dmso_8.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt"

## proteomics 24
samples=norm_proteins[,c(1:3,13:15)]
cl=c(rep(1,3),rep(0,3))
# proteomics 24 hours
out_dir="/home/bianca/Desktop/rp_prot_mkc_vs_dmso_24.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_24_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_24_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_24_DEregulatedBio.txt"

## proteomics 48
samples=norm_proteins[,c(4:6,16:18)]
cl=c(rep(1,3),rep(0,3))
# proteomics 48 hours
out_dir="/home/bianca/Desktop/rp_prot_mkc_vs_dmso_48.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_48_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_48_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_48_DEregulatedBio.txt"

## proteomics 72
samples=norm_proteins[,c(7:9,19:21)]
cl=c(rep(1,3),rep(0,3))
# proteomics 72 hours
out_dir="/home/bianca/Desktop/rp_prot_mkc_vs_dmso_72.txt"
upsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_72_upregulatedBio.txt"
downsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_72_downregulatedBio.txt"
allsave= "/home/bianca/Desktop/rp_prot_mkc_vs_dmso_72_DEregulatedBio.txt"

### # RP parameters:
is_logged <- as.logical("FALSE")
is_logged <- as.logical("TRUE")
logbase <- 2
logbase <- 10
is_normalized <- as.logical("FALSE")
is_normalized <- as.logical("TRUE")
data_size=as.logical("FALSE")

#interconnected
out_dir="/home/bianca/Desktop/interconnected_hypoxia_mkc_4_NEW.txt"
upsave= "/home/bianca/Desktop/interconnected_hypoxia_mkc_4_upregulatedBio_NEW.txt"
downsave= "/home/bianca/Desktop/interconnected_hypoxia_mkc_4_downregulatedBio_NEW.txt"
allsave= "/home/bianca/Desktop/interconnected_hypoxia_mkc_4_DEregulatedBio_NEW.txt"

# hypoxia mkc
out_dir="/home/bianca/Desktop/hypoxia_mkc_4_NEW.txt"
upsave= "/home/bianca/Desktop/hypoxia_mkc_4_upregulatedBio_NEW.txt"
downsave= "/home/bianca/Desktop/hypoxia_mkc_4_downregulatedBio_NEW.txt"
allsave= "/home/bianca/Desktop/hypoxia_mkc_4_DEregulatedBio_NEW.txt"

#hypoxia dmso
out_dir="/home/bianca/Desktop/hypoxia_dmso_4_NEW.txt"
upsave= "/home/bianca/Desktop/hypoxia_dmso_4_upregulatedBio_NEW.txt"
downsave= "/home/bianca/Desktop/hypoxia_dmso_4_downregulatedBio_NEW.txt"
allsave= "/home/bianca/Desktop/hypoxia_dmso_4_DEregulatedBio_NEW.txt"




######### limma edger blah blah
y1<-y
adjusted<- adjusted - min(adjusted)
y <- DGEList(counts=adjusted, samples=sample_info , genes=rownames(adjusted))

y <- estimateDisp(y, mod1)

# mean-variance plot. raw variances of the counts (grey dots), the variances using the tagwise

meanVarPlot <- plotMeanVar( y , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

### DEG
fit <- glmFit(y, model.matrix(~0 + y$samples$treatmentL))
D8vsM8 <- glmLRT(fit, contrast = c(0,-1,0,1))
D24vM24 <- glmLRT(fit, contrast = c(-1,0,1,0))

#message("data_size argument: in case the user is interested in specific num.gene as output")

# if (data_size){
#   num.gene=as.numeric(args[11])  
# }else{
#   
#   method <- as.character(args[11])
#   #message("method can be 'pval', 'pfp'")
#   cutoff <- as.numeric(args[12])
# }



if (data_size){
  num.gene=as.numeric(2000)  
}else{
  
  method <- as.character("pval")
  #message("method can be 'pval', 'pfp'")
  cutoff <- as.numeric(0.05)
}


#FUNCTIONS:

#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  
  return(norm)
}


norm_samples = function.norm(samples)

#RP function:  ###  BECAUSE RP DOES 0 VS 1 0 = condition!!!
function.rp <-function(norm_samples){
  
  RP=RankProd::RankProducts(norm_samples, cl=cl, logged =is_logged, na.rm = TRUE, 
                            gene.names = rownames(norm_samples),
                            plot = FALSE, rand = 123, calculateProduct = TRUE, MinNumOfValidPairs = NA,
                            RandomPairs = NA, huge = FALSE)
  
  if (data_size){
    
    RP_top=RankProd::topGene(RP,gene.names=rownames(samples), num.gene = num.gene, logged = is_logged, logbase = logbase) 
    # logbase=2 is the default
    
  }else{
    
    RP_top=RankProd::topGene(RP,gene.names=rownames(norm_samples), method = method, cutoff = cutoff, logged = is_logged, logbase = logbase) 
    
    if (is.null(rownames(RP_top$Table1)) || is.null(rownames(RP_top$Table2)))
      RP_top1=matrix(rep("NAN",5),nrow=1,ncol=5)
    RP_top2=matrix(rep("NAN",5),nrow=1,ncol=5)
    #print(RP_top2)
    print("Rank Products did not identify differentially expressed features, either up or downregulated")
    
    
    if (!is.null(rownames(RP_top$Table1)))
      RP_top1=as.data.frame(RP_top$Table1) 
    RP_top1[,3]<- (log2(RP_top1[,3]))  ###  BECAUSE RP DOES 0 VS 1 0 = condition!!!
    colnames(RP_top2)<-colnames(RP_top1)
    if (!is.null(rownames(RP_top$Table2))){
      RP_top2=as.data.frame(RP_top$Table2)
      RP_top2[,3]<- (log2(RP_top2[,3])) ### BECAUSE RP DOES 0 VS 1 
      colnames(RP_top1)<-colnames(RP_top2)
    }
    
    
    
    RP_TOP<-rbind(RP_top1, RP_top2)
    
    return(list(RP_TOP, RP))
  }
}



# function for ensembl annotation (preprocessing when rownames are of this type: ENSG00000109320.13)  
# # this will not be used yet
# fun.symbol=function(table){
#   z=c()
#   for (i in rownames(table)){
#     a=strsplit(i,"[[:punct:]]")[[1]][1]
#     z=append(z, a)
#   }
#   return(z)
# }


### analysis

if (!is_normalized){
  # normalize data
  print("Normalizing Data")
  
  function.norm <- function(counts){
    norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
    rownames(norm)=rownames(counts)
    colnames(norm)=colnames(counts)
    
    return(norm)
  }
  
  norm_samples = function.norm(samples)
  
  
  # RankProducts run
  print("RankProducts run")
  RP_TOP=function.rp(norm_samples)[1]
  
  print="\n\ntopGene function run,\nT1: downregulated genes in condition vs control\nT2: upregulated genes in condition vs control\n\n"
  cat(print)
  
  # Saving outputs 
  print("Saving outputs")
  
  #  
  RP_TOP<-as.data.frame(RP_TOP)
  #  rownames(RP_TOP)<-RP_TOP[,1]
  #a<-fun.symbol(RP_TOP)
  #print(head(a))
  #RP_TOP<-cbind(a,as.data.frame(RP_TOP))
  #RP_TOP<-RP_TOP[,c(1:6)]
  write.table(RP_TOP, out_dir,quote=FALSE, sep="\t", dec=".")
  
  print("Worflow Finished")
  
}else{
  
  print("RankProducts run with normalized data")
  
  RP_TOP=function.rp(samples)[1]
  
  print="\n\ntopGene function run,\nT1: downregulated genes in condition vs control\nT2: upregulated genes in condition vs control\n\n"
  cat(print)
  
  # Saving outputs 
  print("Saving outputs")
  RP_TOP<-as.data.frame(RP_TOP)
  #  rownames(RP_TOP)<-RP_TOP[,1]
  #a<-fun.symbol(RP_TOP)
  #print(head(a))
  #RP_TOP<-cbind(a,as.data.frame(RP_TOP))
  #RP_TOP<-RP_TOP[,c(1:6)]
  write.table(RP_TOP, out_dir, quote=FALSE, sep="\t", dec=".")
  
  print("Worflow Finished")
  print("Saving results")
  upregulated<-RP_TOP[RP_TOP$FC..class1.class2. > 0 & RP_TOP$P.value < 0.05,]
  downregulated<-RP_TOP[RP_TOP$FC..class1.class2. < 0 & RP_TOP$P.value < 0.05,]
  write.table(upregulated[,c(3,5)], upsave, sep="\t", quote=F, row.names = rownames(upregulated))
  write.table(RP_TOP[,c(3,5)], allsave, sep="\t", quote=F, row.names = rownames(RP_TOP))
  write.table(downregulated[,c(3,5)], downsave, sep="\t", quote=F, row.names = rownames(downregulated))
  }




# RP_Y<-rownames(RP_TOP)
# length(which(RP_TOP$FC..class1.class2. < -0.5))
# length(which(RP_TOP$FC..class1.class2. > 0.5))