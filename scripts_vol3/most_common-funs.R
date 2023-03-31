library("sva")
library(tidyverse)
library(ggplot2)
library(preprocessCore)
library(RankProd)
library(biomaRt)
library("edgeR")
library("limma")
library("RColorBrewer")
library("stringr")
library("reshape2")


# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp_12 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#882255","#661100","#999933","#44AA99")

### functions batch correction

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

#sampleInfo = DATA FRAME samples treatment replicate
comsva<- function(samples,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  adjustedc <- ComBat(samples, batch=sampleInfo$replicate)
  svaob<-adjustedc - min(adjustedc)
  svobj <- svaseq(svaob, mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(cleaned_count)
}

svacom<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(as.matrix(svaob), mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv[c(1,2)])) 
  adjustedc <- ComBat(cleaned_count, batch=sampleInfo$replicate)
  return(adjustedc)
}

svAll<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(svaob, mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(cleaned_count)
}

com<- function(samples,sampleInfo,replicate){
  adjustedc <- ComBat(samples, batch=sampleInfo$replicate)
  return(adjustedc)
}


### functions pca

PreparePlotpca<-function(counts,sample_info,var1,var2,values){
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
    ggplot(aes(x =PC1, y = PC2, col = var1, shape = var2)) +
    scale_colour_manual(values=values) +
    geom_point(size=5) +
    geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle("PCA of treatment and replicate") +
    theme_bw()
  p
}


## without manual colour
PreparePlotpca2<-function(counts,sample_info,var1,var2,values){
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
    ggplot(aes(x =PC1, y = PC2, col = var1, shape = var2)) +
    scale_fill_brewer() +
    geom_point(size=5) +
    geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle("PCA of treatment and replicate") +
    theme_bw()
  p
}

## without manual colour and no label
PreparePlotpca3<-function(counts,sample_info,var1,var2,values){
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
    ggplot(aes(x =PC1, y = PC2, col = var1, shape = var2)) +
    scale_fill_brewer() +
    geom_point(size=5) +
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle("PCA of treatment and replicate") +
    theme_bw()
  p
}


### title, colour manual, label adjusted (geom_text)
PreparePlotpca6<-function(counts,sample_info,var1,var2,values,title){
  pca_matrix <- counts %>% 
    as.matrix() %>% 
    t()
  sample_pca <- prcomp(pca_matrix)
  tibpca<-as_tibble(pca_matrix, rownames = "sample")
  pc_eigenvalues <- sample_pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>%
    mutate(pct = variance/sum(variance)*100) %>% 
    mutate(pct_cum = cumsum(pct))
  pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct), fill="deeppink4") +
    geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")
  pc_scores <- sample_pca$x
  pc_scores <- pc_scores %>% 
    as_tibble(rownames = "sample")
  #wow<-list(sample_pca,pc_eigenvalues)
  p<-sample_pca$x %>% 
    as_tibble(rownames = "sample") %>% 
    full_join(sample_info, by = "sample") %>% 
    ggplot(aes(x =PC1, y = PC2, col = var1, shape = var2)) +
    scale_fill_brewer() +
    geom_point(size=5) + 
    geom_text(aes(label = sample),vjust = 0.2,nudge_y = 0.2) + 
    labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
    ggtitle(title) +
    theme_bw()
  p
}

############## RANKPROD ############
### # RP parameters:
is_logged <- as.logical("FALSE")
is_logged <- as.logical("TRUE")
logbase <- 2
logbase <- 10
is_normalized <- as.logical("FALSE")
is_normalized <- as.logical("TRUE")
data_size=as.logical("FALSE")

### SAVE OUTPUTS
out_dir="/home/ben/Desktop/time_series_proteomics_rankprod.txt"
upsave= "/home/ben/Desktop/time_series_proteomics_rankprod.txt"
downsave= "/home/ben/Desktop/time_series_proteomics_rankprod.txt"
allsave= "/home/ben/Desktop/time_series_proteomics_rankprod.txt"


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
  RP_TOP<-as.data.frame(RP_TOP)
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
  write.table(RP_TOP, out_dir, quote=FALSE, sep="\t", dec=".")
  print("Worflow Finished")
  print("Saving results")
  upregulated<-RP_TOP[RP_TOP$FC..class1.class2. > 0 & RP_TOP$P.value < 0.05,]
  downregulated<-RP_TOP[RP_TOP$FC..class1.class2. < 0 & RP_TOP$P.value < 0.05,]
  write.table(upregulated[,c(3,5)], upsave, sep="\t", quote=F, row.names = rownames(upregulated))
  write.table(RP_TOP[,c(3,5)], allsave, sep="\t", quote=F, row.names = rownames(RP_TOP))
  write.table(downregulated[,c(3,5)], downsave, sep="\t", quote=F, row.names = rownames(downregulated))
}


