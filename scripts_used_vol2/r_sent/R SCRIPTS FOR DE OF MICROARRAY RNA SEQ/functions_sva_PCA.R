library("sva")
library(tidyverse)
library(ggplot2)


# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp_12 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#882255","#661100","#999933","#44AA99")

######## !!!!!!!!!!!!!  CAREFULL TREATMENT IS TIMEPOINTS*DRUG, CONDITION=DRUG ALONE!!!!!!!!!!!!!!!!!!! ########################### 
########################## YOU CAN EITHER DO TREATMENT OR CONDITION * TIME WHEN YOU TRY TO CORRECT ##########
######### BETTER ALL TOGETHER TO BE SURE OF YOUR DESIGN MATRIX ##############################################


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


### isoform GENE COUNTS ALL

samples
sample_info


# runs
svacomre<-svacom(as.matrix(samples),sample_info,treatment*time,replicate)
comsvare<-comsva(as.matrix(samples),sample_info,treatment*time,replicate)
comre<-com(samples,sample_info,replicate)
savre<-svAll(as.matrix(samples),sample_info,treatment*time,replicate)

## pca plots
#RAW
PreparePlotpca(samples,sample_info,sample_info$sampleId,sample_info$replicate,values=Set3)

# CLEAN
PreparePlotpca(svacomre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(comsvare,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(comre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(savre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)


# ### only 8 hours
# samples8<-samples[,c(1,2,3,4,5,6)]
# sample_info8<-sample_info[c(1:6),]
# 
# # runs
# svacomre8<-svacom(samples8,sample_info8,treatment,replicate)
# comsvare8<-comsva(samples8,sample_info8,treatment,replicate)
# comre8<-com(samples8,sample_inf8o,replicate)
# savre8<-svAll(samples8,sample_info8,treatment,replicate)
# 
# PreparePlotpca(svacomre8,sample_info8,sample_info8$sampleId,sample_info8$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(comsvare8,sample_info8,sample_info8$sampleId,sample_info8$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(comre8,sample_info8,sample_info8$sampleId,sample_info8$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(svare8,sample_info8,sample_info8$sampleId,sample_info8$replicate,values=c("deepskyblue4","deeppink4"))
# 
# ### obly 24 hours
# samples24<-samples[,c(7:12)]
# sample_info24<-sample_info[c[7:12],]
# 
# # runs
# svacomre24<-svacom(samples24,sample_info24,treatment,replicate)
# comsvare24<-comsva(samples24,sample_info24,treatment,replicate)
# comre24<-com(samples24,sample_info24,replicate)
# savre24<-svAll(samples24,sample_info24,treatment,replicate)
# 
# 
# PreparePlotpca(svacomre24,sample_info24,sample_info24$sampleId,sample_info24$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(comsvare24,sample_info24,sample_info24$sampleId,sample_info24$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(comre24,sample_info24,sample_info24$sampleId,sample_info24$replicate,values=c("deepskyblue4","deeppink4"))
# PreparePlotpca(svare24,sample_info24,sample_info24$sampleId,sample_info24$replicate,values=c("deepskyblue4","deeppink4"))

## abundance
samples<-as.data.frame(stringtie$abundance)
rownames(samples)<-samples[,1]
samples<-samples[,-1]
samples<-data.frame(samples[,c(4,5,6,10,11,12,1,2,3,7,8,9)])

# sample_info<-as.data.frame(sampleInfo)
# sample_info$replicate<-c(paste("R",1:3,sep=""))
# sample_info$treatment<-c(rep("DMSO",6),rep("MKC",6))
# sample_info$time<-rep(c(rep(8,3),rep(24,3)),2)
# colnames(sample_info)<-c("sampleId", "sample","replicate","treatment","time")
# sample_info<-sample_info[c(1,2,3,7,8,9,4,5,6,10,11,12),]


### MICROARRAY AND HYPOXIA AT 4 HOURS
# runs
samples<-y$counts
sample_info<-samples_info
sample_info$replicate<-as.factor(sample_info$replicate)
sample_info$treatment<-as.factor(sample_info$treatment)
sample_info$sampleId<-c(rep("DMSO_NORMAL",3),rep("MKC_NORMAL",3),rep("DMSO_HYPOXIA",3),rep("MKC_HYPOXIA",3))

svacomre<-svacom(as.matrix(samples),sample_info,treatment,replicate)
comsvare<-comsva(as.matrix(samples),sample_info,treatment*condition,replicate)
comre<-com(samples,sample_info,replicate)
savre<-svAll(as.matrix(samples),sample_info,treatment*condition,replicate)

## pca plots
#RAW
PreparePlotpca(samples,sample_info,sample_info$sample,sample_info$replicate,values=cbp2)
PreparePlotpca2(samples,sample_info,sample_info$sampleId,sample_info$replicate,values=Set3)
# CLEAN
PreparePlotpca(svacomre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(comsvare,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(comre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)
PreparePlotpca(savre,sample_info,sample_info$sampleId,sample_info$replicate,values=cbp2)

## hif
#RAW
samples<-as.data.frame(expHif)
sample_info<-samples_info
PreparePlotpca(samples,sample_info,sample_info$sample,sample_info$replicate,values=cbp2)
PreparePlotpca2(samples,sample_info,sample_info$sampleId,sample_info$replicate,values=Set3)



#############3 proteomics
samples<-proteins[,c(1:3,10:15,22:24)]
sample_info=sample_info_proteins[c(1:3,10:15,22:24),]
samples=proteins
sample_info=sample_info_proteins
colnames(sample_info)=c("name","condition","hours","replicate","sample")
sample_info$joint=paste(sample_info$condition,sample_info$hours, sep="_")
PreparePlotpca2(samples,sample_info,sample_info$joint,sample_info$replicate,values=Set3)
PreparePlotpca3(samples,sample_info,sample_info$joint,sample_info$replicate,values=Set3)

svacomre<-svacom(as.matrix(samples),sample_info,condition*time,replicate)
comsvare<-comsva(as.matrix(samples),sample_info,condition*hours,replicate)
comre<-com(samples,sample_info,replicate)
savre<-svAll(as.matrix(samples),sample_info,condition*hours,replicate)
