library("sva")
library(tidyverse)
library(ggplot2)
library(preprocessCore)
library(biomaRt)
library("edgeR")
library("limma")
library("RColorBrewer")
library("stringr")
library(optparse)

############################################################ BACK END #######################################################################
#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  return(norm)
}

########################################################################### batch effect corrections

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
  comsvare <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(list(comsvare,svobj))
}


svacom<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(as.matrix(svaob), mod1, mod0) 
  cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv[c(1,2)])) 
  svacomre <- ComBat(cleaned_count, batch=sampleInfo$replicate)
  return(list(svacomre,svobj))
}

svAll<- function(svaob,sampleInfo,treatment,replicate){
  mod0 <- model.matrix(~1, as.data.frame(sampleInfo))
  mod1 <- model.matrix(~0 + treatment, as.data.frame(sampleInfo))
  svobj <- svaseq(svaob, mod1, mod0) 
  svAllre <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 
  return(list(svAllre,svobj))
}

com<- function(samples,sampleInfo,replicate){
  comre <- ComBat(samples, batch=sampleInfo$replicate)
  return(comre)
}
### functions pca
# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# colourblindpalette
cbp_12 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#882255","#661100","#999933","#44AA99")

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
################# prepare DGElist and preprocess
preproc<-function(counts,si,design,rpnorm){
  ### make object, filter and normalize
  counts <- mutate_all(counts, function(x) as.numeric(as.character(x)))
  y <- DGEList(counts=counts, samples=si, genes=rownames(counts))
  # REMOVE LOW EXPRESSED GENES
  keep.exprs <- filterByExpr(y,design=design, min.count=3)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  ####### normalization
  if(rpnorm){
    y$counts<-function.norm(y$counts)
  }
  else{
    y <- calcNormFactors(y, method = "upperquartile")
  }
  return(y)
}


difglmQ<-function(y,design,contr.matrix,pval){
  y <- estimateDisp(y, design, robust = T)
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmQLFit(y, design)
  resall <- glmQLFTest(fit, contrast=contr.matrix)
  fres<-topTags(resall, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = pval)
  return(list(resall,fres))
}

difglmQA<-function(y,design,contr.matrix,pval){
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  resall <- glmQLFTest(fit, contrast=contr.matrix)
  fres<-topTags(resall, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = pval)
  return(list(resall,fres))
}

difglmF<-function(y,design,contr.matrix,pval){
  y <- estimateDisp(y, design, robust = T)
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  resall <- glmLRT(fit, contrast=contr.matrix)
  fres<-topTags(resall, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = pval)
  return(list(resall,fres))
}

difclasslm<-function(counts,design,contr.matrix,lfc,coef){
  v<-voom(counts,design)
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  tfit <- treat(vfit, lfc=lfc)
  res <- topTreat(tfit, coef=coef, n=Inf)
  return(res)
}


#####################################################################################################################################################




counts<-read.table("/home/ben/Desktop/pacli/all_together_pac.txt")
si<-read.table("/home/ben/Desktop/pacli/all_together_pac_samplesInfo.txt")
si2<-si


design <- model.matrix(~0 + groups, si)
colnames(design)<-make.names(colnames(design))
contr.matrix8 <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  levels = colnames(design))
contr.matrix24 <- makeContrasts(
  pacdmso8 = groupspac24-groupsdmso24, 
  levels = colnames(design))
contr.matrix <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  pacdmso24 = groupspac24 - groupsdmso24,
  levels = colnames(design))

y<-preproc(counts,si,design,FALSE)

pac8<-difglm(y,design,contr.matrix8, 0.05)
pac24<-difglm(y,design,contr.matrix24, 0.05)

################
#batch effect correction and visualization
colnames(y$counts)=si$sample
PreparePlotpca(y$counts, si, si$treatment, as.factor(si$replicate), cbp_12)
comsvaresres<-comsva(y$counts, si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svacomrere<-svacom(y$counts, si , treatment, replicate)
PreparePlotpca(svacomrere[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svaresres<-svAll(y$counts, si , treatment, replicate)
PreparePlotpca(svaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
comresres<-com(y$counts, si , replicate)
PreparePlotpca(comresres, si, si$treatment, as.factor(si$replicate), cbp_12)


################## 8
si<-si[c(1:6),]

PreparePlotpca(y$counts[,c(1:6)], si, si$treatment, as.factor(si$replicate), cbp_12)
comsvaresres<-comsva(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svacomrere<-svacom(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(svacomrere[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svaresres<-svAll(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(svaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
comresres<-com(y$counts[,c(1:6)], si , replicate)
PreparePlotpca(comresres, si, si$treatment, as.factor(si$replicate), cbp_12)

################## 24
si<-si2[c(7:12),]
PreparePlotpca(y$counts[,c(7:12)], si, si$treatment, as.factor(si$replicate), cbp_12)
comsvaresres<-comsva(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svacomrere<-svacom(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(svacomrere[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svaresres<-svAll(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(svaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
comresres<-com(y$counts[,c(7:12)], si , replicate)
PreparePlotpca(comresres, si, si$treatment, as.factor(si$replicate), cbp_12)
#####################################

colnames(counts)<-si$sample
y<-preproc(counts,si,design,FALSE)
comsvaresres<-comsva(y$counts, si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)


design1 <- cbind(design,comsvaresres[[2]][[1]])
colnames(design1)<-make.names(colnames(design1))

contr.matrix8 <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  levels = colnames(design1))
contr.matrix24 <- makeContrasts(
  pacdmso8 = groupspac24-groupsdmso24, 
  levels = colnames(design1))
contr.matrix <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  pacdmso24 = groupspac24 - groupsdmso24,
  levels = colnames(design1))

y<-preproc(counts,si,design1,FALSE)
pac8a<-difglmQ(y,design1,contr.matrix8, 0.05)
pac24a<-difglmQ(y,design1,contr.matrix24, 0.05)
paca<-difglmQ(y,design1,contr.matrix, 0.05)

pac8af<-difglmF(y,design1,contr.matrix8, 0.05)
pac24af<-difglmF(y,design1,contr.matrix24, 0.05)
pacaf<-difglmF(y,design1,contr.matrix, 0.05)




table(rowSums(y$counts==0)==12)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
y$samples$norm.factors





