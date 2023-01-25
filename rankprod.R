#!/usr/bin/env Rscript

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

NOT YET INSTALLED 
library(pd.hugene.2.0.st)
#library(affy)
#library(oligo)
library(affycoretools)




### inputs microarray
gse <- getGEO("GSE99766", GSEMatrix = TRUE)
#columns<-pData(gse[[1]])
#fData(gse[[1]])
exp<-exprs(gse[[1]])

#### dmso 4 h with hypoxia 10, 11, 12 columns
expIn<-exp[,1:6]
expInrev<-exp[,c(4,5,6,1,2,3)]
#colnames(expIn)<-c(rep(0,3), rep(1,3))
annotation(gse) <- "hugene20sttranscriptcluster.db" # check if necessary

# optional for alternative ways


gse1<- annotateEset(gse, hugene20sttranscriptcluster.db)

samples_vector_file<-data.frame(c(rep(0,3), rep(1,3)))
colnames(samples_vector_file)<- "condition"
samples_vector_file$replicate<-c(rep(c("1","2","3"), 2))
#samples_vector_file$replicate<-c(rep(c("a","b","c"), 2))
samples_vector_file$sample <-c(paste("DMSO", 1:3, sep=""), paste("MKC",1:3, sep=""))
samples_vector_file$treatment <-c(rep("DMSO", 3), rep("MKC",3))
colnames(expIn)<-samples_vector_file$sample

#annotation of genes
httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))

listAttributes(mart)
vaslues=rownames(samples)

homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")

ncbi_to_gene<-getBM(attributes = c("entrezgene_id","hgnc_symbol"), 
                    filters="entrezgene_id", values=rownames(expIn), mart=mart)


genes <- getBM(attributes = c("affy_huex_1_0_st_v2", "hgnc_symbol"), filters = "affy_huex_1_0_st_v2", values = rownames(samples), mart = mart)




# create DGEobject 
y <- DGEList(counts=expIn, samples=samples_vector_file, genes=rownames(expIn))
names(y)


# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]


# COUNT PER MILLION, LOG2 COUNTS. for plotting 
lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
y$counts <- cpm(y) 


## Normalisation by the method of trimmed mean of M-values (TMM)

nsamples <- ncol(expIn)
col <- brewer.pal(nsamples, "Paired")

wee <- log2(y$counts)
boxplot(wee, las=2, col=col, main="")
title(main="Log2 Raw data",ylab="Log-cpm")

boxplot(lcpm_pre, las=2, col=col, main="")
title(main="Log2 CPM data",ylab="Log-cpm")

y <- calcNormFactors(y)
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")



### EXPLORATORY ANALYSIS - VISUALIZATION

### NULL MODEL ONLY AGJUSTMENT VARIABLES (not the ones that we are interested in)
mod0 <- model.matrix(~0 + repeats, y$samples)
#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + repeats + group, data=y$samples)


#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod1 <- model.matrix(~1, y$samples)

#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod0 <- model.matrix(~0 + group, data=y$samples)


#### svobj = sva(edata,mod,mod0,n.sv=n.sv)
svobj <- svaseq(cpm(y), mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.
log_cleaned_count <- log2(cleaned_count)

#
boxplot(log_cleaned_count, las=2, col=col, main="")
title(main="SVA normalised data",ylab="Log-cpm")

### PCA to inspect batch correction
# no sva
pca <- prcomp(t(na.omit(y$counts)), center = T, scale. = T)

summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= samples_vector_file$pca1)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

ggplot(PCAi , aes(PC1, PC2, col= samples_vector_file$replicate)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by replicate") +
  theme_bw() 

ggplot(PCAi , aes(PC2, PC3, col= samples_vector_file$pca1)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC3", round(eigs[3] / sum(eigs), digits = 2)), x = paste("PC2", round(eigs[2] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

ggplot(PCAi , aes(PC2, PC3, col= samples_vector_file$replicate)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC3", round(eigs[3] / sum(eigs), digits = 2)), x = paste("PC2", round(eigs[2] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by replicate") +
  theme_bw() 


ggplot(PCAi , aes(PC3, PC4, col= samples_vector_file$pca1)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC4", round(eigs[4] / sum(eigs), digits = 2)), x = paste("PC3", round(eigs[3] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

ggplot(PCAi , aes(PC3, PC4, col= samples_vector_file$replicate)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC4", round(eigs[4] / sum(eigs), digits = 2)), x = paste("PC3", round(eigs[3] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by replicate") +
  theme_bw() 

# after sva
pca <- prcomp(t(na.omit(log_cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)

variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 


head(log_cleaned_count)

################# RNA seq experiment

mrna<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/mRNA/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/Data/raw_counts.txt", sep="\t", row.names=1)
raw_counts <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/raw_counts.txt"), sep="", quote = "FALSE") 

sample_info_RNAseq <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/sample_info_RNAseq.txt"), sep="") 

eset <-raw_counts[,7:18] # first 6 columns is annotation

# 
# gse <- getGEO("GSE176454", GSEMatrix = TRUE)
# columns<-pData(gse[[1]])
# fData(gse[[1]])
# exp<-exprs(gse[[1]])

### my own preparation
# mrna1<-mrna[,(7:18)]
# mrnakeep<-mrna1[,c(1,5,4,9,8,10,3,2,6,12,11,7)]
# rmrna<-mrnakeep
# colnames(rmrna)<-c(paste("DMSO_8",1:3,sep="_"), paste("MKC_8", 1:3,sep="_"), paste("DMSO_24",1:3,sep="_"), paste("MKC_24",1:3,sep="_"))
# samples_vector_file<-data.frame(c(rep(c(rep(0,3), rep(1,3)),2)))
# colnames(samples_vector_file)<- "condition"
# samples_vector_file$replicate<-rep(c(as.numeric(rep(c("1","2","3"), 2))),2)
# #maybe change pca a little to show time points as well ?
# samples_vector_file$pca <-rep(c(paste("DMSO", 1:3, sep=""), paste("MKC",1:3, sep="")),2)



#Annotate genes
listEnsembl(GRCh=38)
listEnsembl(version=91)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
todo_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), mart=ensembl)

genetable <- data.frame(gene.id=rownames(eset))
genes <- merge(genetable, todo_genes, by.x= "gene.id", by.y="ensembl_gene_id", all.x=TRUE)

dup <- genes$gene.id[duplicated(genes$gene.id)]
mat <- match(genetable$gene.id, genes$gene.id)
genes <- genes[mat,]

# create DGEobject
y <- DGEList(counts=eset, samples=sample_info_RNAseq, genes=rownames(eset))
names(y)


# create DGEobject
y <- DGEList(counts=eset, samples=sample_info_RNAseq, genes=genes)
names(y)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]



# COUNT PER MILLION, LOG2 COUNTS. for plotting
cpm <- cpm(y) 
lcpm_pre <- cpm(y, log=TRUE)

## Normalisation by the method of trimmed mean of M-values (TMM)

nsamples <- ncol(eset)
col <- brewer.pal(nsamples, "Paired")

wee <- log2(y$counts)
boxplot(wee, las=2, col=col, main="")
title(main="Log2 Raw data",ylab="Log-cpm")

boxplot(lcpm_pre, las=2, col=col, main="")
title(main="Log2 CPM data",ylab="Log-cpm")

y <- calcNormFactors(y)
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

######### BATCH CORRECTON AND EXPLORATORY ANALYSIS
### PCA to inspect batch correction
# no sva
pca <- prcomp(t(na.omit(y$counts)), center = T, scale. = T)

summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# after sva
pca <- prcomp(t(na.omit(log_cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)

variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 


head(log_cleaned_count)

# BATCH EFFECT CORRECTION WITH SVA

### NULL MODEL ONLY AGJUSTMENT VARIABLES (not the ones that we are interested in)
mod0 <- model.matrix(~0 + repeats, y$samples)
#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + repeats + group, data=y$samples)


#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod1 <- model.matrix(~1, y$samples)

#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod0 <- model.matrix(~0 + group, data=y$samples)


#### svobj = sva(edata,mod,mod0,n.sv=n.sv)
svobj <- svaseq(cpm(y), mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.
log_cleaned_count <- log2(cleaned_count)

#
boxplot(log_cleaned_count, las=2, col=col, main="")
title(main="SVA normalised data",ylab="Log-cpm")


### COMBAT 

adjusted <- ComBat_seq(lcpm, batch=y$samples$repeats, group=NULL)



### PCA to inspect batch correction
# no sva
pca <- prcomp(t(na.omit(lcpm_pre)), center = T, scale. = T)

summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# after sva
pca <- prcomp(t(na.omit(log_cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)

variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 


head(log_cleaned_count)


PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$repeats)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by repeat") +
  theme_bw() 



#### pca after combat

pca <- prcomp(t(na.omit(adjusted)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)

variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")



PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$treatment))  + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment after combat") +
  theme_bw() 




##### DIFFERENTIAL EXPRESSION WITH RANKPROD
#my rna proscess
samples=rmrna

#microarray
samples=(expIn)
cl<-samples_vector_file$condition

# rna seq sample info
sample_info_RNAseq$rankprod<-c(rep(1,6), rep(0,6))

#rna seq mrna 8 h
samples=log_cleaned_count[,c(1,4,5,8,9,10)]
cl<-sample_info_RNAseq$rankprod[c(1,4,5,8,9,10)]

#rna seq mrna 24 h
samples=log_cleaned_count[,c(2,3,6,7,11,12)]
cl<-sample_info_RNAseq$rankprod[c(2,3,6,7,11,12)]

### parameters
is_logged <- as.logical("FALSE")
is_logged <- as.logical("TRUE")
logbase <- 2
is_normalized <- as.logical("FALSE")
is_normalized <- as.logical("TRUE")
out_dir="/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/mRNA/rankprod_results_24_rev.txt"

# RP parameters:
data_size=as.logical("FALSE")
#message("data_size argument: in case the user is interested in specific num.gene as output")

if (data_size){
  num.gene=as.numeric(args[11])  
}else{
  
  method <- as.character(args[11])
  #message("method can be 'pval', 'pfp'")
  cutoff <- as.numeric(args[12])
}



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

#RP function: 
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
    RP_top1[,3]<- log2(RP_top1[,3])
    colnames(RP_top2)<-colnames(RP_top1)
    if (!is.null(rownames(RP_top$Table2))){
      RP_top2=as.data.frame(RP_top$Table2)
      RP_top2[,3]<- log2(RP_top2[,3])
      colnames(RP_top1)<-colnames(RP_top2)
    }
    
    
    
    RP_TOP<-rbind(RP_top1, RP_top2)
    
    return(list(RP_TOP, RP))
  }
}



# function for ensembl annotation (preprocessing when rownames are of this type: ENSG00000109320.13)  
# this will not be used yet
fun.symbol=function(table){
  z=c()
  for (i in rownames(table)){
    a=strsplit(i,"[[:punct:]]")[[1]][1]
    z=append(z, a)
  }
  return(z)
}


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
  
}
