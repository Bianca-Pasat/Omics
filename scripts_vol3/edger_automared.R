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
  mod1 <- model.matrix(~0 + sampleInfo$treatment, as.data.frame(sampleInfo))
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
  y <- DGEList(counts=counts, samples=si)
  # REMOVE LOW EXPRESSED GENES
  keep.exprs <- filterByExpr(y,design=design)
  y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  ####### normalization
  if(rpnorm){
    y$counts<-function.norm(y$counts)
  }
  else if(!rpnorm){
    y <- calcNormFactors(y, method = "upperquartile")
  }
  return(y)
}

boxbo<-function(counts,sample_info,hours,treatment, title){
  raw_cts_long <- counts %>% 
    pivot_longer(1:length(colnames(counts)), names_to = "sample", values_to = "expression")
  
  # Join with sample information table
  raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
  
  raw_cts_long %>%
    ggplot(aes(as.character(hours), log2(expression),fill = treatment)) + 
    geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
    facet_grid(cols = vars(replicate)) +
    labs (y = "log2 expression", x = "time points") +
    ggtitle(title) +
    theme_bw()
  
}

### TODO put variables
densbo<-function(counts,sample_info,title){
  raw_cts_long <-counts %>% 
    pivot_longer(1:length(colnames(counts)), names_to = "sample", values_to = "expression")
  # Join with sample information table
  raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
  raw_cts_long %>%
         ggplot( aes(x=log2(expression))) +
         geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
         ggtitle(title) +
         theme_bw()

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


heatmapBianca<-function(paclimeta, genesign){
  testx<-paclimeta[,c(1,2)]
  rownames(testx) <- rownames(paclimeta)
  colnames(testx) <- colnames(paclimeta[,c(1,2)])
  testx<-testx[genesign,]
  test<-na.omit(testx)
  col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
  fontsize=0.5
  signa<-as.data.frame(paclimeta[,c(3:length(colnames(paclimeta)))])
  col<-HeatmapAnnotation(df=paclimeta[,c(3:length(colnames(paclimeta)))])
  row<-HeatmapAnnotation(signature=signa ,which="row")
  hm<-Heatmap(as.matrix(test), name="expression",
              col=col_fun,
              cluster_rows=T,
              cluster_columns = T,
              #row_names_side = "left",
              show_row_names = F,
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
              ##bottom_annotation = col,
              heatmap_width = unit(15, "cm"),
              width = NULL,
              heatmap_height = unit(15, "cm"))
  draw(hm)
  
}

testx<-paclimeta[,c(1,2)]
rownames(testx) <- rownames(paclimeta)
colnames(testx) <- colnames(paclimeta[,c(1,2)])
testx<-testx[genesign,]
test<-na.omit(test)
test<-testx
col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5
signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(df=paclimeta[,c(3:length(colnames(paclimeta)))])
row<-HeatmapAnnotation(signature=signa[,1] ,which="row")



#####################################################################################################################################################




counts<-read.table("/home/ben/Desktop/pacli/all_together_pac.txt")
si<-read.table("/home/ben/Desktop/pacli/all_together_pac_samplesInfo.txt")
colnames(counts)<-si$sample
si2<-si


design <- model.matrix(~0 + groups, si)
colnames(design)<-make.names(colnames(design))
boxbo(counts,si,hours,treatment,"Boxplot of raw counts before fitlering")
y<-preproc(counts,si,design,TRUE)
boxbo(y$counts,si,hours,treatment,"Boxplot of raw counts after normalization")



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

pac8<-difglmF(y,design,contr.matrix8, 0.05)
pac24<-difglmF(y,design,contr.matrix24, 0.05)
pacNoclean<-difglmF(y,design,contr.matrix, 0.05)
################
#batch effect correction and visualization
siAll<-si
siAll$treat<-si2$treatment
siAll$treatment<-si2$groups
colnames(y$counts)=siAll$sample
PreparePlotpca(y$counts, siAll, siAll$treatment, as.factor(siAll$replicate), cbp_12)
comsvaresres<-comsva(y$counts, siAll , groups, replicate)
comsvaresresT<-comsva(y$counts, siAll , treatment, replicate)
PreparePlotpca(comsvaresresT[[1]], siAll, siAll$treatment, as.factor(siAll$replicate), cbp_12)
PreparePlotpca(comsvaresresT[[1]], siAll, siAll$groups, as.factor(siAll$replicate), cbp_12)
PreparePlotpca(comsvaresres[[1]], siAll, siAll$groups, as.factor(siAll$replicate), cbp_12)
svacomrere<-svacom(y$counts, siAll , treatment, replicate)
PreparePlotpca(svacomrere[[1]], siAll, siAll$treatment, as.factor(siAll$replicate), cbp_12)
svaresres<-svAll(y$counts, siAll , treatment, replicate)
PreparePlotpca(svaresres[[1]], siAll, siAll$treatment, as.factor(siAll$replicate), cbp_12)
comresres<-com(y$counts, siAll , replicate)
PreparePlotpca(comresres, siAll, siAll$treatment, as.factor(siAll$replicate), cbp_12)


################## 8
si<-si2[c(1:6),]

PreparePlotpca(y$counts[,c(1:6)], si, si$treatment, as.factor(si$replicate), cbp_12)
comsvaresres8<-comsva(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svacomrere8<-svacom(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(svacomrere8[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svaresres8<-svAll(y$counts[,c(1:6)], si , treatment, replicate)
PreparePlotpca(svaresres8[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
comresres8<-com(y$counts[,c(1:6)], si , replicate)
PreparePlotpca(comresres8, si, si$treatment, as.factor(si$replicate), cbp_12)

################## 24
si<-si2[c(7:12),]
PreparePlotpca(y$counts[,c(7:12)], si, si$treatment, as.factor(si$replicate), cbp_12)
comsvaresres24<-comsva(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(comsvaresres24[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svacomrere24<-svacom(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(svacomrere24[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
svaresres24<-svAll(y$counts[,c(7:12)], si , treatment, replicate)
PreparePlotpca(svaresres24[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)
comresres24<-com(y$counts[,c(7:12)], si , replicate)
PreparePlotpca(comresres24, si, si$treatment, as.factor(si$replicate), cbp_12)
#####################################

colnames(counts)<-si$sample
y<-preproc(counts,si,design,FALSE)
comsvaresres<-comsva(y$counts, si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, as.factor(si$replicate), cbp_12)

design <- model.matrix(~0 + groups, si)
design1 <- cbind(design,comsvaresres[[2]][[1]])
colnames(design1)<-make.names(colnames(design1))


#### pacli with cleaned on all
y$counts <- svaob
design1<-cbind(design, svobj$sv)
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

#pac alone with sv
counts<-counts[,c(1:6)]
y<-preproc(counts,si,design1,FALSE)
pac8a<-difglmQ(y,design1,contr.matrix8, 0.05)



table(rowSums(y$counts==0)==12)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
y$samples$norm.factors


atf4<-read.table("/home/ben/Desktop/ATF_sign.txt", sep="\t")
paclimeta<-as.data.frame(paca[[2]][,c(1,2)])
paclimeta$hif<-ifelse(rownames(paclimeta) %in% hif$V1, "HIF","no")
paclimeta$senescence<-ifelse(rownames(paclimeta) %in% sen$V1, "Senesence","no")
paclimeta$xbp1<-ifelse(rownames(paclimeta) %in% xbp1$V1, "XBP1","no")
paclimeta$hours<-ifelse(rownames(paclimeta) %in% rownames(pac8a[[2]]), "8 hours","24 hours")
paclimeta$atf4<-ifelse(rownames(paclimeta) %in% xbp1$V1, "ATF4","no")

write.table(paclimeta,"/home/ben/Desktop/pacli/palimeta_all_de.txt",sep="\t", quote=F)

sign_xbp1=c("ASS1", "C3", "CCL20", "COL4A6", "CXCL2", "CXCL5", "CXCL8", "IFI44L", "IL1B", "IL6", "KCNN2", "MMP1", "MMP12", "MMP3", "PLA2G4A", "PPP4R4", "SERPINB2", "TFPI2",
            "ZNF804A")
priori8<-read.table("/home/ben/Desktop/pacli/8HOURS_PRIORI_of_cleanedAll.tsv", sep="\t")
priori24<-read.table("/home/ben/Desktop/pacli/24_hours_24priori_of_cleanedAll.tsv", sep="\t")
genesign<-c(priori8$V2, priori24$V2)
genesign<-genesign[-1]
paclimeta$priori<-ifelse(rownames(paclimeta) %in% genesign , "Prioritized","no")
vec<-c(which(paclimeta$atf4!="no"), which(paclimeta$senescence!="no"),which(paclimeta$hif!="no"),which(paclimeta$priori!="no"))

newPac<-paclimeta[vec,]
newPac<-newPac[-c(22,25,29,32,36,38,39),]
newPac$names<-rownames(newPac)
n8<-newPac[which(newPac$hours=="8 hours"),]
n8<-n8[,-c(2)]
n24<-newPac[which(newPac$hours=="24 hours"),]
n24<-n24[,-c(1)]

newlong8<-n8 %>% pivot_longer(c(2,3,4,6,7),  values_to="features")
newlong24<-n24 %>% pivot_longer(c(2,3,4,6,7),  values_to="features")

newlong24$features<-gsub("no", NA, newlong24$features)
newlong8$features<-gsub("no", NA, newlong8$features)
newl8<-na.omit(newlong8)
ggplot(newl24,aes(names,logFC.pacdmso24,color = features)) + 
       geom_boxplot() +  
       theme_bw()
ggplot(newl8,aes(names,logFC.pacdmso8,color = features)) + 
  geom_boxplot() +  
  theme_bw()

heatmapBianca(paclimeta, rownames(paclimeta))
heatmapBianca(paclimeta, genesign)


features