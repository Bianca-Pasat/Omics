
#BiocManager::install(c("plotrix","STATegRa","MASS","gridExtra","r.jive"))
#BiocManager::install(c("RegularizedSCA","RpESCA","RSpectra"))
# library(devtools)
# install_github('YipengUva/RpESCA')
library(preprocessCore)
library("plotrix")
library("STATegRa")
library("MASS")
library("gridExtra")
library("r.jive")
library("RegularizedSCA")
library("RpESCA")
library("RSpectra")
library(plyr)
#read rna, mirna, proteomics


#Joint omics exploration CAREFULL ON THE PC S YOURS BIGHT BE DIFFERENT!!!!!!!!!!!!!1
#Common samples
# Load require packages
library(VennDiagram)
library(plotrix)
# Venn Diagram
m <- list(rna=colnames(rna_gbm), 
          mirna=colnames(mirna_gbm), 
          met=colnames(logit_met_gbm))

vennPlot <- venn.diagram(m, NULL, fill=c('red','green','blue'), 
                         cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2)
grid.newpage()
grid.draw(vennPlot)

#OmicsPCA exploration

# Load require packages
library(STATegRa)
library(MASS)
library(gridExtra)

#Data scaling

#Calculate Forbenious normalization
frobenius_rna <- norm(as.matrix(log2_rna_gbm), type="F")
frobenius_mirna <- norm(as.matrix(mirna_gbm), type="F")
frobenius_met <- norm(as.matrix(log_proteins), type="F")

#Mean-centring and division by Frobenius normalization factor
rna_gbm_norm <- t(scale(t(rna_gbm), scale=FALSE))/frobenius_rna
mirna_gbm_norm <- t(scale(t(mirna_gbm), scale=FALSE))/frobenius_mirna
met_gbm_norm <- t(scale(t(logit_met_gbm), scale=FALSE))/frobenius_met

##Mean-centring without Frobenius normalization
rna_gbm_norm2 <- t(scale(t(rna_gbm), scale=FALSE))
mirna_gbm_norm2 <- t(scale(t(mirna_gbm), scale=FALSE))
met_gbm_norm2 <- t(scale(t(logit_met_gbm), scale=FALSE))

#Common samples between mRNA and miRNA
ids_rna_mirna <- intersect(colnames(rna_gbm), colnames(mirna_gbm)) 
length(ids_rna_mirna)


#Common samples between mRNA and Methylation
ids_rna_met <- intersect(colnames(rna_gbm), colnames(logit_met_gbm))
length(ids_rna_met)
##model selection

# JIVE
# Load require packages
library(r.jive)

#Using r.jive

####################################################
## mRNA + miRNA (already centered and normalized) ##
####################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_mirna],
             miRNA=mirna_gbm_norm[,ids_rna_mirna])
Data<-mda
rjive_rna_mirna <- jive(Data, scale = FALSE, center = FALSE) 


##########################################################
## mRNA + methylation (already centered and normalized) ##
##########################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_met],
             met=met_gbm_norm[,ids_rna_met])
rjive_rna_met <- jive(Data, scale = FALSE, center = FALSE)

#PCA-GCA
# Load require packages
library(RegularizedSCA)

#Using pca-gca

####################################################
## mRNA + miRNA (already centered and normalized) ##
####################################################
data <- cbind(t(rna_gbm_norm[,ids_rna_mirna]), 
              t(mirna_gbm_norm[,ids_rna_mirna]))

data<-cbind(t(rna_gbm_norm), t(mirna_gbm_norm))
# Number of variables
Jk <- c(nrow(rna_gbm_norm), nrow(mirna_gbm_norm))

# Select components
pca_gca(data, Jk, cor_min = .7)

##########################################################
## mRNA + methylation (already centered and normalized) ##
##########################################################
data <- cbind(t(rna_gbm_norm[,ids_rna_met]), 
              t(met_gbm_norm[,ids_rna_met]))

# Number of variables
Jk <- c(nrow(rna_gbm_norm), nrow(met_gbm_norm))

# Select components
pca_gca(data, Jk, cor_min = .7)

#pESCA
# Load require packages
library(RpESCA)
library(RSpectra)

#Using RpESCA

################################################
## mRNA + miRNA (already centered)            ##
################################################

dataSets <- list(rna=t(rna_gbm_norm2[,ids_rna_mirna]), 
                 mirna=t(mirna_gbm_norm2[,ids_rna_mirna]))

dataSets<-list(rna=t(rna_gbm_norm), mirna=t(mirna_gbm_norm))
dataTypes <- 'GG'

# used parameters
opts <- list()
opts$tol_obj <- 1E-4
opts$quiet <- 1

# save the estimated alphas, selected number of PCs and cvErrors
alphas <- rep(NA,2)
R_selected_list <- as.list(1:2)
cvErrors_list <- as.list(1:2)

# alpha estimation for each data set CAREFULL! Rs is dimention of data............)
for (i in 1:length(dataSets)){
  # index ith data set
  X <- dataSets[[i]]
  
  # alpha estimation procedure
  alpha_est <- alpha_estimation(X,K=2,Rs =1:6,opts=opts)
  alphas[i] <- alpha_est$alphas_mean
  R_selected_list[[i]] <- alpha_est$R_CV
  cvErrors_list[[i]] <- alpha_est$cvErrors
}
names(alphas) <- paste0('alpha', 1:2)

## Setting the parameters for the pESCA model

# concave function and its hyper-parameter
fun_concave <- 'gdp'; gamma <- 1;
penalty = 'L2' # concave L2 norm penalty

# Parameters of a pESCA with concave L2norm penalty model
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-6 # stopping criteria
opts$maxit   <- 500
opts$alphas  <- alphas
#opts$R # components used Default 0.5*min(I,J) - Minimum number of variables.
opts$thr_path <- 0 # generaint thresholding path or not
opts$quiet <- 1

# Model selection pESCA conave L2 norm penalty
nTries <- 15
lambdas_CV <- log10_seq(from=1, to=500, length.out=nTries)

result_CV <- pESCA_CV(dataSets, dataTypes,
                      lambdas_CV, 
                      penalty=penalty, 
                      fun_concave=fun_concave, 
                      opts=opts)


# select the model with minimum CV error
index_min_cv <- which.min(result_CV$cvErrors_mat[,1])

# cvErrors during the model selection process
result_CV$cvErrors_mat

# selected value of lambda
lambdas_CV[index_min_cv]

# fit the final model
lambdas_opt <- rep(lambdas_CV[index_min_cv],length(dataSets))
opts_opt <- result_CV$inits[[index_min_cv]]
opts_opt$tol_obj <- 1E-8 # using high precision model


pESCA_L2 <- pESCA(dataSets = dataSets,
                  dataTypes = dataTypes,
                  lambdas = lambdas_opt,
                  penalty = penalty,
                  fun_concave= fun_concave,
                  opts = opts_opt)


mu <- pESCA_L2$mu
A <- pESCA_L2$A
B <- pESCA_L2$B
S <- pESCA_L2$S

# estimated variation explained ratios for each data set or the full data set
pESCA_L2$varExpTotals

# estimated variation explained ratios of each PC for each data set or the full data set
pESCA_L2$varExpPCs

#Cut-off 1% 
sel1 <- which(pESCA_L2$varExpPCs[1,]>1) #14
sel2 <- which(pESCA_L2$varExpPCs[2,]>1) #23
intersect(sel1,sel2)

#Cut-off 5%
sel1 <- which(pESCA_L2$varExpPCs[1,]>10) #3
sel2 <- which(pESCA_L2$varExpPCs[2,]>10) #2
intersect(sel1,sel2)


################################################
## mRNA + Methylation (already centered)      ##
################################################

dataSets <- list(rna=t(rna_gbm_norm2[,ids_rna_met]), 
                 met=t(met_gbm_norm2[,ids_rna_met]))

dataTypes <- 'GG'

# used parameters
opts <- list()
opts$tol_obj <- 1E-4
opts$quiet <- 1

# save the estimated alphas, selected number of PCs and cvErrors
alphas <- rep(NA,2)
R_selected_list <- as.list(1:2)
cvErrors_list <- as.list(1:2)

# alpha estimation for each data set
for (i in 1:length(dataSets)){
  # index ith data set
  X <- dataSets[[i]]
  
  # alpha estimation procedure
  alpha_est <- alpha_estimation(X,K=3,Rs = 5:20,opts=opts)
  alphas[i] <- alpha_est$alphas_mean
  R_selected_list[[i]] <- alpha_est$R_CV
  cvErrors_list[[i]] <- alpha_est$cvErrors
}
names(alphas) <- paste0('alpha', 1:2)

## Setting the parameters for the pESCA model

# concave function and its hyper-parameter
fun_concave <- 'gdp'; gamma <- 1;
penalty = 'L2' # concave L2 norm penalty

# Parameters of a pESCA with concave L2norm penalty model
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-6 # stopping criteria
opts$maxit   <- 500
opts$alphas  <- alphas
#opts$R  # components used Default 0.5*min(I,J) - Minimum number of variables.
opts$thr_path <- 0 # generaint thresholding path or not
opts$quiet <- 1

# Model selection pESCA conave L2 norm penalty
nTries <- 15
lambdas_CV <- log10_seq(from=1, to=500, length.out=nTries)

result_CV <- pESCA_CV(dataSets, dataTypes,
                      lambdas_CV, 
                      penalty=penalty, 
                      fun_concave=fun_concave, 
                      opts=opts)


# select the model with minimum CV error
index_min_cv <- which.min(result_CV$cvErrors_mat[,1])

# cvErrors during the model selection process
result_CV$cvErrors_mat

# selected value of lambda
lambdas_CV[index_min_cv]

# fit the final model
lambdas_opt <- rep(lambdas_CV[index_min_cv],length(dataSets))
opts_opt <- result_CV$inits[[index_min_cv]]
opts_opt$tol_obj <- 1E-8 # using high precision model


pESCA_L2 <- pESCA(dataSets = dataSets,
                  dataTypes = dataTypes,
                  lambdas = lambdas_opt,
                  penalty = penalty,
                  fun_concave= fun_concave,
                  opts = opts_opt)

mu <- pESCA_L2$mu
A <- pESCA_L2$A
B <- pESCA_L2$B
S <- pESCA_L2$S

# estimated variation explained ratios for each data set or the full data set
pESCA_L2$varExpTotals

# estimated variation explained ratios of each PC for each data set or the full data set
pESCA_L2$varExpPCs

#Cut-off 1%
sel1 <- which(pESCA_L2$varExpPCs[1,]>1) #3
sel2 <- which(pESCA_L2$varExpPCs[2,]>1) #16
intersect(sel1,sel2)
# 2 common components
# 1 dist mRNA + 14 methylation
# Total: 17 components

#Cut-off 5%
sel1 <- which(pESCA_L2$varExpPCs[1,]>5) #3
sel2 <- which(pESCA_L2$varExpPCs[2,]>5) #3
intersect(sel1,sel2)

#Sub-space recovery
#Based on the number of components identified (pESCA) we can recover the subspace.
#It is interesting to take a look on the common and distinctive components identified through this method.

################################################
## mRNA + miRNA (already centered)            ##
################################################

# Explore the common and distinctive components and those variation explained
# Common components
cc <- intersect(sel1,sel2)
pESCA_L2$varExpPCs[,cc]


# Distinctive components
dist_rna <- sel1[!sel1%in%cc] #rna
pESCA_L2$varExpPCs[,dist_rna]

dist_mirna <- sel2[!sel2%in%cc] #mirna
pESCA_L2$varExpPCs[,dist_mirna]

# Keep de common and distinctive components
common_pc <- pESCA_L2$A[,cc]
colnames(common_pc) <- colnames(pESCA_L2$varExpPCs[,cc])

dist_rna_pc <- pESCA_L2$A[,dist_rna]
colnames(dist_rna_pc) <- names(dist_rna)

dist_mirna_pc <- pESCA_L2$A[,dist_mirna]
colnames(dist_mirna_pc) <- names(dist_mirna)


# metadata for common samples
clin <- clinical_gbm[ids_rna_mirna,]

# PCA plot
require(ggplot2)

df <- as.data.frame(common_pc)
df$group <- clin$GeneExp_Subtype
df$batch <- clin$batch_number

#plot of first common components by gene expression subtype
p <- ggplot(df,aes_string("PC2","PC3",color="group"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by Gene Expression Subtype") + theme(legend.position="bottom")


#plot of first common components by batch
p <- ggplot(df,aes_string("PC2","PC3",color="batch"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by batch") + theme(legend.position="bottom")


#Let’s explore how relate these common and distinctive components with the clinical variables.

df <- data.frame(common_pc,dist_rna_pc,dist_mirna_pc)
colnames(df) <- c(paste("common_",colnames(common_pc),sep=""),
                  paste("dist_rna_",colnames(dist_rna_pc),sep=""),
                  paste("dist_mirna_",colnames(dist_mirna_pc),sep=""))

pcaCorrelation(clinical=clin, clin_var, df, maxlogPvalue =10, maxPC=ncol(df), pc.names=TRUE)



################################################
## mRNA + methylation (already centered)      ##
################################################

# Explore the common and distinctive components and those variation explained
# Common components
cc <- intersect(sel1,sel2)
pESCA_L2$varExpPCs[,cc]

# Distinctive components
dist_rna <- sel1[!sel1%in%cc] #rna
pESCA_L2$varExpPCs[,dist_rna]
##        X_1        X_2     X_full 
## 13.0520328  0.0000000  0.7036662
dist_met <- sel2[!sel2%in%cc] #met
pESCA_L2$varExpPCs[,dist_met]

# Keep de common and distinctive components
common_pc <- pESCA_L2$A[,cc]
colnames(common_pc) <- colnames(pESCA_L2$varExpPCs[,cc])

dist_rna_pc <- data.frame(pESCA_L2$A[,dist_rna])
colnames(dist_rna_pc) <- names(dist_rna)

dist_met_pc <- pESCA_L2$A[,dist_met]
colnames(dist_met_pc) <- names(dist_met)


# metadata for common samples
clin <- clinical_gbm[ids_rna_met,]

# PCA plot
require(ggplot2)

df <- as.data.frame(common_pc)
df$group <- clin$GeneExp_Subtype
df$batch <- clin$batch_number

#plot of first common components by gene expression subtype
p <- ggplot(df,aes_string("PC2","PC7",color="group"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by Gene Expression Subtype") + theme(legend.position="bottom")


#plot of first common components by batch
p <- ggplot(df,aes_string("PC2","PC7",color="batch"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by batch") + theme(legend.position="bottom")



df <- data.frame(common_pc,dist_rna_pc,dist_met_pc)
colnames(df) <- c(paste("common_",colnames(common_pc),sep=""),
                  paste("dist_rna_",colnames(dist_rna_pc),sep=""),
                  paste("dist_met_",colnames(dist_met_pc),sep=""))

pcaCorrelation(clinical=clin, clin_var, df, maxlogPvalue =10, maxPC=ncol(df), pc.names=TRUE)



#Batch correction

# Load require package
library(sva)


summary(as.factor(as.character(clinical_gbm[,'batch_number'])))

clinical_gbm[which(clinical_gbm$batch_number%in%c('26.41.0','62.42.0')),c(1,4)]

rm_samples_batch <- c('TCGA-27-1836-01', 'TCGA-32-2495-01')

rna_gbm_corrected <- rna_gbm[,!colnames(rna_gbm)%in%rm_samples_batch] 


mirna_gbm_corrected <- mirna_gbm[,!colnames(mirna_gbm)%in%rm_samples_batch] 

#Expression data (rna_gbm_corrected)

mod <- model.matrix(~gender+GeneExp_Subtype+days_to_birth, 
                    data=clinical_gbm[colnames(rna_gbm_corrected),])

batch <- as.factor(as.character(clinical_gbm[colnames(rna_gbm_corrected),
                                             'batch_number']))

rna_gbm_corrected <- ComBat(as.matrix(rna_gbm_corrected),
                            batch=batch,
                            mod=mod, 
                            prior.plot=TRUE)


#Exploring the PCA from corrected data we can see that the batch effect have been removed.

pca_rna <- prcomp(t(rna_gbm_corrected))
var_prcomp <- pca_rna$sdev^2

# Plot PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), pc=c(1:length(var_prcomp)))
ggplot(pcvar[1:15,], aes(x = pc)) + geom_line(aes(y=var)) + 
  labs(x='Principal Component', y='Explained Variance') + 
  geom_point(aes(y=var)) + 
  ggtitle('mRNA: Principal Components variance')


# Plot PCs scores
df <- data.frame(PC1=pca_rna$x[,1], PC2=pca_rna$x[,2],
                 clinical_gbm[colnames(rna_gbm_corrected),
                              c("GeneExp_Subtype","gender","batch_number")])

ggplot(df) +
  geom_point(aes(x=PC1, y=PC2,color=factor(GeneExp_Subtype)),size=5,shape=20) +
  guides(color=guide_legend("GeneExp_Subtype"), 
         fill=guide_legend("GeneExp_Subtype")) +
  labs(x='PC1 (12.8%)', y='PC2 (8.4%)') +
  ggtitle("Principal Components from batch-adjusted mRNA data") +
  theme(legend.position="bottom")


ggplot(df) +
  geom_point(aes(x=PC1, y=PC2, color=factor(batch_number)), size=5, shape=20) +
  guides(color=guide_legend("batch_number"),
         fill=guide_legend("batch_number")) +
  labs(x='PC1 (12.8%)', y='PC2 (8.4%)') +
  ggtitle("Principal Components from batch-adjusted mRNA data") +
  theme(legend.position="bottom")


#Now, the association analysis shows how the batch effect is clearly removed from our data 
# and the main variability is related with gene expression subtype.

pcaCorrelation(clinical=clinical_gbm[colnames(rna_gbm_corrected),],
               clin_var=clin_var, pca_data=pca_rna$x, maxlogPvalue = 10)


#miRNA data (mirna_gbm_corrected)
mod <- model.matrix(~gender+GeneExp_Subtype+days_to_birth, 
                    data=clinical_gbm[colnames(mirna_gbm_corrected),])

batch <- as.factor(as.character(clinical_gbm[colnames(mirna_gbm_corrected),
                                             'batch_number']))

mirna_gbm_corrected <- ComBat(as.matrix(mirna_gbm_corrected),
                              batch=batch,
                              mod=mod, 
                              prior.plot=TRUE)


#Exploring the PCA from corrected data we can see that the batch effect have been removed. 
# However, in the case of miRNA, the results are equivalent to before batch adjustment, 
# as no batch-effect was detected in the initial data set.

#miRNA
pca_mirna <- prcomp(t(mirna_gbm_corrected))
var_prcomp <- pca_mirna$sdev^2


# Plot PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), pc=c(1:length(var_prcomp)))
ggplot(pcvar[1:15,], aes(x = pc)) + geom_line(aes(y=var)) + 
  labs(x='Principal Component', y='Explained Variance') + 
  geom_point(aes(y=var)) + 
  ggtitle('miRNA: Principal Components variance')


# Plot PCs scores
df <- data.frame(PC1=pca_mirna$x[,1], PC2=pca_mirna$x[,2],
                 clinical_gbm[colnames(mirna_gbm_corrected),
                              c("GeneExp_Subtype","gender","batch_number")])

ggplot(df) +
  geom_point(aes(x=PC1, y=PC2, color=factor(GeneExp_Subtype)),size=5,shape=20) +
  guides(color=guide_legend("GeneExp_Subtype"),
         fill=guide_legend("GeneExp_Subtype")) +
  labs(x='PC1 (19.5%)', y='PC2 (9.5%)') +
  ggtitle("Principal Components from batch-adjusted miRNA data") +
  theme(legend.position="bottom")


ggplot(df) +
  geom_point(aes(x=PC1, y=PC2, color=factor(batch_number)), size=5, shape=20) +
  guides(color=guide_legend("batch_number"),
         fill=guide_legend("batch_number")) +
  labs(x='PC1 (19.5%)', y='PC2 (9.5%)') +
  ggtitle("Principal Components from batch-adjusted miRNA data") +
  theme(legend.position="bottom")


pcaCorrelation(clinical=clinical_gbm[colnames(mirna_gbm_corrected),],
               clin_var=clin_var, pca_data=pca_mirna$x, maxlogPvalue = 10)


#Methylation data (logit_met_gbm)


mod <- model.matrix(~gender+GeneExp_Subtype+days_to_birth, 
                    data=clinical_gbm[colnames(logit_met_gbm),])

batch <- as.factor(as.character(clinical_gbm[colnames(logit_met_gbm),
                                             'batch_number']))

logit_met_gbm_corrected <- ComBat(logit_met_gbm,
                                  batch=batch, 
                                  mod=mod, 
                                  prior.plot=TRUE, 
                                  par.prior=FALSE, 
                                  BPPARAM=MulticoreParam())


pca_met <- prcomp(t(logit_met_gbm_corrected))
var_prcomp <- pca_met$sdev^2

# Plot PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), pc=c(1:length(var_prcomp)))
ggplot(pcvar[1:15,], aes(x = pc)) + geom_line(aes(y=var)) + 
  labs(x='Principal Component', y='Explained Variance') + 
  geom_point(aes(y=var)) + 
  ggtitle('Methylation: Principal Components variance')


# Plot PCs scores
df <- data.frame(PC1=pca_met$x[,1], PC2=pca_met$x[,2],
                 clinical_gbm[colnames(logit_met_gbm_corrected),
                              c("GeneExp_Subtype","gender","batch_number")])

ggplot(df) +
  geom_point(aes(x=PC1, y=PC2,color=factor(GeneExp_Subtype)),size=5,shape=20) +
  guides(color=guide_legend("GeneExp_Subtype"),
         fill=guide_legend("GeneExp_Subtype")) +
  labs(x='PC1 (14.8%)', y='PC2 (11.7%)') +
  ggtitle("Principal Components from batch-adjusted Methylation data") +
  theme(legend.position="bottom")


ggplot(df) +
  geom_point(aes(x=PC1, y=PC2, color=factor(batch_number)), size=5, shape=20) +
  guides(color=guide_legend("batch_number"),
         fill=guide_legend("batch_number")) +
  labs(x='PC1 (14.8%)', y='PC2 (11.7%)') +
  ggtitle("Principal Components from batch-adjusted Methylation data") +
  theme(legend.position="bottom")


#Now, the association analysis shows how the batch effect is clearly removed from our data and the main variability is related with gene expression subtype and survival outcome.

pcaCorrelation(clinical=clinical_gbm[colnames(logit_met_gbm_corrected),],
               clin_var=clin_var, pca_data=pca_met$x, maxlogPvalue = 10)


#Characterization of the data: Joint exploration
#Common samples

# Load require packages
library(VennDiagram)
library(plotrix)
# Venn Diagram
m <- list(rna=colnames(rna_gbm_corrected), 
          mirna=colnames(mirna_gbm_corrected), 
          met=colnames(logit_met_gbm_corrected))

vennPlot <- venn.diagram(m, NULL, fill=c('red','green','blue'), 
                         cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2)
grid.newpage()
grid.draw(vennPlot)



#OmicsPCA exploration


# Load require packages
library(STATegRa)
library(MASS)
library(gridExtra)

#Calculate Forbenious normalization
frobenius_rna <- norm(as.matrix(rna_gbm_corrected), type="F")
frobenius_mirna <- norm(as.matrix(mirna_gbm_corrected), type="F")
frobenius_met <- norm(as.matrix(logit_met_gbm_corrected), type="F")

#Mean-centring and division by Frobenius normalization factor
rna_gbm_norm <- t(scale(t(rna_gbm_corrected), scale=FALSE))/frobenius_rna
mirna_gbm_norm <- t(scale(t(mirna_gbm_corrected), scale=FALSE))/frobenius_mirna
met_gbm_norm <- t(scale(t(logit_met_gbm_corrected), scale=FALSE))/frobenius_met

##Mean-centring without Frobenius normalization
rna_gbm_norm2 <- t(scale(t(rna_gbm_corrected), scale=FALSE))
mirna_gbm_norm2 <- t(scale(t(mirna_gbm_corrected), scale=FALSE))
met_gbm_norm2 <- t(scale(t(logit_met_gbm_corrected), scale=FALSE))


#Common samples between mRNA and miRNA
ids_rna_mirna <- intersect(colnames(rna_gbm_corrected), colnames(mirna_gbm_corrected)) 
length(ids_rna_mirna)

#Common samples between mRNA and Methylation
ids_rna_met <- intersect(colnames(rna_gbm_corrected), colnames(logit_met_gbm_corrected))
length(ids_rna_met)

#Model selection
#Now, we are ready to explore if there are share components between data-types.

#The first step is to determine the number of share/common and individual components between omics.

#How to determine these numbers is still an open question nowadays, so the user can explore different methods like jive, pca-gca and pESCA, between others. Following, we show how to perform the model selection based on these three methods.

#JIVE
# Load require packages
library(r.jive)
#Using r.jive

################################################
## mRNA + miRNA (already centered)            ##
################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_mirna],
             miRNA=mirna_gbm_norm[,ids_rna_mirna])
rjive_rna_mirna <- jive(Data, scale = FALSE, center = FALSE) 


######################################################
## mRNA + methylation (already centered)            ##
######################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_met],
             met=met_gbm_norm[,ids_rna_met])
rjive_rna_met <- jive(Data, scale = FALSE, center = FALSE)

#PCA-GCA
# Load require packages
library(RegularizedSCA)

#Using pca-gca

################################################
## mRNA + miRNA (already centered and scaled) ##
################################################
data <- cbind(t(rna_gbm_norm[,ids_rna_mirna]), 
              t(mirna_gbm_norm[,ids_rna_mirna]))

# Number of variables
Jk <- c(nrow(rna_gbm_norm), nrow(mirna_gbm_norm))

# Select components
pca_gca(data, Jk, cor_min = .7)

######################################################
## mRNA + methylation (already centered and scaled) ##
######################################################
data <- cbind(t(rna_gbm_norm[,ids_rna_met]), 
              t(met_gbm_norm[,ids_rna_met]))

# Number of variables
Jk <- c(nrow(rna_gbm_norm), nrow(met_gbm_norm))

# Select components
pca_gca(data, Jk, cor_min = .7)
# 2.6.2.2.3pESCA
# Load require packages
library(RpESCA)
library(RSpectra)

#Using RpESCA

################################################
## mRNA + miRNA (already centered)            ##
################################################

dataSets <- list(rna=t(rna_gbm_norm2[,ids_rna_mirna]), 
                 mirna=t(mirna_gbm_norm2[,ids_rna_mirna]))

dataTypes <- 'GG'

# used parameters
opts <- list()
opts$tol_obj <- 1E-4
opts$quiet <- 1

# save the estimated alphas, selected number of PCs and cvErrors
alphas <- rep(NA,2)
R_selected_list <- as.list(1:2)
cvErrors_list <- as.list(1:2)

# alpha estimation for each data set
for (i in 1:length(dataSets)){
  # index ith data set
  X <- dataSets[[i]]
  
  # alpha estimation procedure
  alpha_est <- alpha_estimation(X,K=3,Rs = 5:20,opts=opts)
  alphas[i] <- alpha_est$alphas_mean
  R_selected_list[[i]] <- alpha_est$R_CV
  cvErrors_list[[i]] <- alpha_est$cvErrors
}
names(alphas) <- paste0('alpha', 1:2)

## Setting the parameters for the pESCA model

# concave function and its hyper-parameter
fun_concave <- 'gdp'; gamma <- 1;
penalty = 'L2' # concave L2 norm penalty

# Parameters of a pESCA with concave L2norm penalty model
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-6 # stopping criteria
opts$maxit   <- 500
opts$alphas  <- alphas
#opts$R  # components used Default 0.5*min(I,J) - Minimum number of variables.
opts$thr_path <- 0 # generaint thresholding path or not
opts$quiet <- 1

# Model selection pESCA conave L2 norm penalty
nTries <- 15
lambdas_CV <- log10_seq(from=1, to=500, length.out=nTries)

result_CV <- pESCA_CV(dataSets, dataTypes,
                      lambdas_CV, 
                      penalty=penalty, 
                      fun_concave=fun_concave, 
                      opts=opts)


# select the model with minimum CV error
index_min_cv <- which.min(result_CV$cvErrors_mat[,1])

# cvErrors during the model selection process
result_CV$cvErrors_mat

# selected value of lambda
lambdas_CV[index_min_cv]

# fit the final model
lambdas_opt <- rep(lambdas_CV[index_min_cv],length(dataSets))
opts_opt <- result_CV$inits[[index_min_cv]]
opts_opt$tol_obj <- 1E-8 # using high precision model


pESCA_L2 <- pESCA(dataSets = dataSets,
                  dataTypes = dataTypes,
                  lambdas = lambdas_opt,
                  penalty = penalty,
                  fun_concave= fun_concave,
                  opts = opts_opt)


mu <- pESCA_L2$mu
A <- pESCA_L2$A
B <- pESCA_L2$B
S <- pESCA_L2$S

# estimated variation explained ratios for each data set or the full data set
pESCA_L2$varExpTotals

# estimated variation explained ratios of each PC for each data set or the full data set
pESCA_L2$varExpPCs

#Cut-off 1% 
sel1 <- which(pESCA_L2$varExpPCs[1,]>1) #13
sel2 <- which(pESCA_L2$varExpPCs[2,]>1) #21
intersect(sel1,sel2)


# Cut-off 5%
sel1 <- which(pESCA_L2$varExpPCs[1,]>5) #5
sel2 <- which(pESCA_L2$varExpPCs[2,]>5) #2
intersect(sel1,sel2)



################################################
## mRNA + Methylation (already centered)      ##
################################################

dataSets <- list(rna=t(rna_gbm_norm2[,ids_rna_met]), 
                 met=t(met_gbm_norm2[,ids_rna_met]))

dataTypes <- 'GG'

# used parameters
opts <- list()
opts$tol_obj <- 1E-4
opts$quiet <- 1

# save the estimated alphas, selected number of PCs and cvErrors
alphas <- rep(NA,2)
R_selected_list <- as.list(1:2)
cvErrors_list <- as.list(1:2)

# alpha estimation for each data set
for (i in 1:length(dataSets)){
  # index ith data set
  X <- dataSets[[i]]
  
  # alpha estimation procedure
  alpha_est <- alpha_estimation(X,K=3,Rs = 5:20,opts=opts)
  alphas[i] <- alpha_est$alphas_mean
  R_selected_list[[i]] <- alpha_est$R_CV
  cvErrors_list[[i]] <- alpha_est$cvErrors
}
names(alphas) <- paste0('alpha', 1:2)

## Setting the parameters for the pESCA model

# concave function and its hyper-parameter
fun_concave <- 'gdp'; gamma <- 1;
penalty = 'L2' # concave L2 norm penalty

# Parameters of a pESCA with concave L2norm penalty model
opts <- list()
opts$gamma <- gamma  # hyper-parameter for the used penalty
opts$rand_start <- 0
opts$tol_obj <- 1e-6 # stopping criteria
opts$maxit   <- 500
opts$alphas  <- alphas
#opts$R # components used Default 0.5*min(I,J) - Minimum number of variables.
opts$thr_path <- 0 # generaint thresholding path or not
opts$quiet <- 1

# Model selection pESCA conave L2 norm penalty
nTries <- 15
lambdas_CV <- log10_seq(from=1, to=500, length.out=nTries)

result_CV <- pESCA_CV(dataSets, dataTypes,
                      lambdas_CV, 
                      penalty=penalty, 
                      fun_concave=fun_concave, 
                      opts=opts)


# select the model with minimum CV error
index_min_cv <- which.min(result_CV$cvErrors_mat[,1])

# cvErrors during the model selection process
result_CV$cvErrors_mat

# selected value of lambda
lambdas_CV[index_min_cv]

# fit the final model
lambdas_opt <- rep(lambdas_CV[index_min_cv],length(dataSets))
opts_opt <- result_CV$inits[[index_min_cv]]
opts_opt$tol_obj <- 1E-8 # using high precision model


pESCA_L2 <- pESCA(dataSets = dataSets,
                  dataTypes = dataTypes,
                  lambdas = lambdas_opt,
                  penalty = penalty,
                  fun_concave= fun_concave,
                  opts = opts_opt)

mu <- pESCA_L2$mu
A <- pESCA_L2$A
B <- pESCA_L2$B
S <- pESCA_L2$S

# estimated variation explained ratios for each data set or the full data set
pESCA_L2$varExpTotals

# estimated variation explained ratios of each PC for each data set or the full data set
pESCA_L2$varExpPCs

#Cut-off 1%
sel1 <- which(pESCA_L2$varExpPCs[1,]>1) #2
sel2 <- which(pESCA_L2$varExpPCs[2,]>1) #9
intersect(sel1,sel2)
# 1 common components
# 1 dist mRNA + 8 methylation
# Total: 10 components

#Cut-off 5%
sel1 <- which(pESCA_L2$varExpPCs[1,]>5) #2
sel2 <- which(pESCA_L2$varExpPCs[2,]>5) #3
intersect(sel1,sel2)

#Sub-space recovery
#Based on the number of components identified (pESCA) we can recover the subspace.

#It is interesting to take a look on the common and distinctive components identified thorugh this method.

################################################
## mRNA + miRNA (already centered)            ##
################################################

# Explore the common and distinctive components and those variation explained
# Common components
cc <- intersect(sel1,sel2)
pESCA_L2$varExpPCs[,cc]

# Distinctive components
dist_rna <- sel1[!sel1%in%cc] #rna
pESCA_L2$varExpPCs[,dist_rna]

dist_mirna <- sel2[!sel2%in%cc] #mirna
pESCA_L2$varExpPCs[,dist_mirna]

# Keep de common and distinctive components
common_pc <- pESCA_L2$A[,cc]
colnames(common_pc) <- colnames(pESCA_L2$varExpPCs[,cc])

dist_rna_pc <- pESCA_L2$A[,dist_rna]
colnames(dist_rna_pc) <- names(dist_rna)

dist_mirna_pc <- pESCA_L2$A[,dist_mirna]
colnames(dist_mirna_pc) <- names(dist_mirna)


# metadata for common samples
clin <- clinical_gbm[ids_rna_mirna,]

# PCA plot
require(ggplot2)

df <- as.data.frame(common_pc)
df$group <- clin$GeneExp_Subtype
df$batch <- clin$batch_number

#plot of first common components by gene expression subtype
p <- ggplot(df,aes_string("PC1","PC3",color="group"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by Gene Expression Subtype") + theme(legend.position="bottom")


#plot of first common components by batch
p <- ggplot(df,aes_string("PC1","PC3",color="batch"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by batch") + theme(legend.position="bottom")


#Let’s explore how relate these common and distinctive components with the clinical variables.

df <- data.frame(common_pc,dist_rna_pc,dist_mirna_pc)
colnames(df) <- c(paste("common_",colnames(common_pc),sep=""),
                  paste("dist_rna_",colnames(dist_rna_pc),sep=""),
                  paste("dist_mirna_",colnames(dist_mirna_pc),sep=""))

pcaCorrelation(clinical=clin, clin_var, df, maxlogPvalue =10, maxPC=ncol(df), pc.names=TRUE)



################################################
## mRNA + methylation (already centered)      ##
################################################

# Explore the common and distinctive components and those variation explained
# Common components
cc <- intersect(sel1,sel2)
pESCA_L2$varExpPCs[,cc]

# Distinctive components
dist_rna <- sel1[!sel1%in%cc] #rna
pESCA_L2$varExpPCs[,dist_rna]

dist_met <- sel2[!sel2%in%cc] #met
pESCA_L2$varExpPCs[,dist_met]

# Keep de common and distinctive components
common_pc <- data.frame(pESCA_L2$A[,cc])
colnames(common_pc) <- paste("PC", cc, sep="")

dist_rna_pc <- data.frame(pESCA_L2$A[,dist_rna])
colnames(dist_rna_pc) <- names(dist_rna)

dist_met_pc <- pESCA_L2$A[,dist_met]
colnames(dist_met_pc) <- names(dist_met)


# metadata for common samples
clin <- clinical_gbm[ids_rna_met,]

# PCA plot
require(ggplot2)

df <- as.data.frame(common_pc)
df$group <- clin$GeneExp_Subtype
df$batch <- clin$batch_number

# As only one common component is found, an auxiliar varaible is defined in order to visualize the first components.
df$auxiliar <- runif(nrow(df),min(df$PC2, na.rm=TRUE),max(df$PC2, na.rm=TRUE))

#plot of first common components by gene expression subtype
p <- ggplot(df,aes_string("PC2","auxiliar",color="group"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by Gene Expression Subtype") + theme(legend.position="bottom")


#plot of first common components by batch
p <- ggplot(df,aes_string("PC2","auxiliar",color="batch"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by batch") + theme(legend.position="bottom")


# Let’s explore how relate these common and distinctive components with the clinical variables.

df <- data.frame(common_pc,dist_rna_pc,dist_met_pc)
colnames(df) <- c(paste("common_",colnames(common_pc),sep=""),
                  paste("dist_rna_",colnames(dist_rna_pc),sep=""),
                  paste("dist_met_",colnames(dist_met_pc),sep=""))

pcaCorrelation(clinical=clin, clin_var, df, maxlogPvalue =10, maxPC=ncol(df), pc.names=TRUE)


#integrative differential analysis by NPC
#Now, we want to identify genes that, according to all modalities considered as a whole (mRNA, miRNA, methylation, etc.), are either deregulated or associated to an outcome of interest.

#From the previous steps (individual omics exploration and joint exploration) we will be considering (1) which co-variates are necessary to include and (2) which analysis do we want to perform.

#Let’s explore the combination of two omics, then NPC will provide two outcomes:
  
#The features that are differentiated in both omics, therefore we are looking into coordinated effects.
#The features that are differentiated in only one of the omics.
#These kind of approach allows different number of samples between omics, so we will explore both escenarios.

#mRNA + miRNA

#Mapping file
#Load the mapping info
#load('hsa-vtm-gene.Rdata')
#the mapping can be obtained using SpidermiR package

#Identify the miRBase version of mapping file and data
library(miRNAmeConverter)

miRNAs <- rownames(mirna_gbm_corrected)
miRNAs_mapping <- names(id)
nc <- MiRNANameConverter()
assessVersion(nc, miRNAs, verbose = FALSE) #version 9.2

assessVersion(nc, miRNAs_mapping, verbose = FALSE) #version 20

#convert miRNA annotation from our miRNA matrix from v.9.2 to v.20
library(anamiR)
new_mirna_matrix <- miR_converter(mirna_gbm_corrected, remove_old = TRUE, original_version=9.2, latest_version = 20)

#restricting mapping info to measured miRNA
id_mapping <- id[intersect(names(id), rownames(new_mirna_matrix))] 

#creating the miRNA association tables
mirna2GeneMap <- stack(id_mapping) #convert in a data.frame. Each line has a gene and one mirna. 
names(mirna2GeneMap) <- c('Gene', 'miRNA')


head(mirna2GeneMap)

#creating the mRNA association tables
expr2GeneMap <- data.frame(measurement = rownames(rna_gbm_corrected), 
                           Gene = rownames(rna_gbm_corrected))


head(expr2GeneMap)

#creating the data mapping
dataMappingExprMirna <- combiningMappings(mappings = list(expr = expr2GeneMap, 
                                                          mirna = mirna2GeneMap),
                                          retainAll = TRUE, reference = 'Gene'); 


head(dataMappingExprMirna)

#We can see that based on our mRNA and miRNA analysis, we have 24,665 pairs in total, it corresponds to 7,814 unique genes (64.9% from the whole mRNA dataset) and 323 unique miRNAs (61.3% from teh whole miRNA dataset).

#Then, we define the relevant clinical variables in each scenario. Our variable of interest in the case of GBM is the survivalOtucome. Moreover, days_to_birth is included as a confounder.

clinicalVariables <- list(expr = c('days_to_birth', 'survivalOutcome'), 
                          met = c('days_to_birth', 'survivalOutcome'),
                          mirna = c('days_to_birth', 'survivalOutcome'));
#Then, we can explore the combination of these two omics using a parametric or non parametric combination (PC or NPC). Moreover, we can perform the analysis on all available samples or just for those that are common between datasets. So, these are the different scenarios that we can have:




#Overlaping samples
#restrict the datasets to the elements of the data mapping and those samples that are common between omics
exprTMP <- rna_gbm_corrected[unique(na.omit(dataMappingExprMirna$expr)), ids_rna_mirna]; #7,814 x 515

mirnaTMP <- new_mirna_matrix[unique(na.omit(dataMappingExprMirna$mirna)), ids_rna_mirna]; #323 x 515

#specifying the data types.
dataTypesExprMirna <- list(ttCoxphContinuous,ttCoxphContinuous)

#preparing the datasets
exprTMP <- createOmicsExpressionSet(Data = as.matrix(exprTMP), 
                                    pData = clinical_gbm[ids_rna_mirna,clinicalVariables[['expr']]])
names(pData(exprTMP)) <- c("age","outcome")

mirnaTMP <- createOmicsExpressionSet(Data = as.matrix(mirnaTMP), 
                                     pData = clinical_gbm[ids_rna_mirna,clinicalVariables[['mirna']]])
names(pData(mirnaTMP)) <- c("age","outcome")

dataInputExprMirna <- list(expr = exprTMP, mirna = mirnaTMP)
#Parametric Combination
set.seed(12345)
omicsPCRes_cPC <- omicsPC(dataInput = dataInputExprMirna,
                          dataMapping = dataMappingExprMirna,
                          dataTypes = dataTypesExprMirna,
                          verbose = TRUE)
#The output of omicsPC provides information regarding p-value and adjusted p-value from both individual and combined analysis.

res <- omicsPCRes_cPC
str(res)

data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$mirna[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$mirna))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$mirna[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$mirna))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot PC
data_plot <- data.frame(pval=unlist(res$pvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$pvaluesPC[,-c(1,2)]),
                                               each=nrow(res$pvaluesPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) + 
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsPC')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$qvaluesPC[,-c(1,2)]),
                                               each=nrow(res$qvaluesPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval,
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


#Now, based on FDR results we can check if there is an increase of statistical power when both omics are jointly analyzed.

#note: here we consider significant an FDR<0.05.

#significant genes
(sig_expr_cPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])

#significant miRNA
(sig_mirna_cPC <- rownames(res$mirna)[which(res$mirna[,"adj.P.Val"]<0.05)])

#Significant genes from mRNA + miRNA parametric combination (Fisher, FDR<0.05)
sig_expr_mirna_cPC <- res$qvaluesPC[which(res$qvaluesPC$Fisher<0.05),1:2] 


#venn diagram
m1 <- list(rna=sig_expr_cPC, rna_fisher=unique(sig_expr_mirna_cPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, 
                          main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(mirna=sig_mirna_cPC, mirna_fisher=unique(sig_expr_mirna_cPC$mirna))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2,
                          main='Significant miRNA identified with individual (mirna) or integration (mirna_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot2)


#When we perform the combination analysis we identify 337 new genes related with GBM survival and 45 miRNAs.
#Here we can see the list of some of the new genes (n=337) identified based on Fisher combination.

setdiff(sig_expr_mirna_cPC$expr,sig_expr_cPC)[1:50]



setdiff(sig_expr_mirna_cPC$mirna,sig_mirna_cPC)

#Non Parametric Combination
set.seed(12345)

# Setting methods to combine pvalues
combMethods <- c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms <- 1000
# Setting number of cores
numCores <- 4
# Setting omicsNPC to print out the steps that it performs.
verbose <- TRUE

omicsNPCRes_cNPC_rna_mirna <- omicsNPC(dataInput = dataInputExprMirna,
                                       dataMapping = dataMappingExprMirna,
                                       dataTypes = dataTypesExprMirna,
                                       combMethods = combMethods, 
                                       numPerms = numPerms,
                                       numCores = numCores,
                                       verbose = verbose)
#The output of omicsNPC provides information regarding p-value and adjusted p-value from both individual and combined analysis.

res <- omicsNPCRes_cNPC_rna_mirna

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$mirna[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$mirna))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$mirna[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$mirna))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot NPC
data_plot <- data.frame(pval=unlist(res$pvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$pvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$pvaluesNPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsNPC')


#FDR plot NPC
data_plot <- data.frame(pval=unlist(res$qvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$qvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$qvaluesNPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsNPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


#Now, based on FDR results we can check if there is an increase of statistical power when both omics are jointly analyzed.


#Significant genes or miRNA from individual omics (FDR<0.05)

#significant genes
(sig_expr_cNPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant miRNA
(sig_mirna_cNPC <- rownames(res$mirna)[which(res$mirna[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + miRNA non parametric combination (Fisher, FDR<0.05)
(sig_expr_mirna_cNPC <- res$qvaluesNPC[which(res$qvaluesNPC[,'Fisher']<0.05),1:2])


#venn diagram
m1 <- list(rna=sig_expr_cNPC, rna_fisher=unique(sig_expr_mirna_cNPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: NPC')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(mirna=sig_mirna_cNPC, mirna_fisher=unique(sig_expr_mirna_cNPC$mirna))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant miRNA identified with individual (mirna) or integration (mirna_fisher)\n omics analysis: NPC')
grid.newpage()
grid.draw(vennPlot2)


#When we perform the combination analysis we identify 23 new genes related with GBM survival and 1 miRNAs.

#Here we can see part of the list of new genes identified based on Fisher combination.

setdiff(sig_expr_mirna_cNPC$expr,sig_expr_cNPC)



setdiff(sig_expr_mirna_cNPC$mirna,sig_mirna_cNPC)

#All samples
#restrict the datasets to the elements of the data mapping
exprTMP <- rna_gbm_corrected[unique(na.omit(dataMappingExprMirna$expr)),]; #7,814 x 523

mirnaTMP <- new_mirna_matrix[unique(na.omit(dataMappingExprMirna$mirna)),]; #323 x 518

#specifying the data types.
dataTypesExprMirna <-  list(ttCoxphContinuous,ttCoxphContinuous)

#preparing the datasets
exprTMP <- createOmicsExpressionSet(Data = as.matrix(exprTMP), 
                                    pData = clinical_gbm[colnames(rna_gbm_corrected),clinicalVariables[['expr']]])
names(pData(exprTMP)) <- c("age","outcome")

mirnaTMP <- createOmicsExpressionSet(Data = as.matrix(mirnaTMP), 
                                     pData = clinical_gbm[colnames(mirna_gbm_corrected),clinicalVariables[['mirna']]])
names(pData(mirnaTMP)) <- c("age","outcome")

dataInputExprMirna <- list(expr = exprTMP, mirna = mirnaTMP)
#Parametric Combination
set.seed(12345)
omicsPCRes_aPC <- omicsPC(dataInput = dataInputExprMirna,
                          dataMapping = dataMappingExprMirna,
                          dataTypes = dataTypesExprMirna,
                          verbose = TRUE)
#The output of omicsPC provides information regarding p-value and adjusted p-value from both individual and combined analysis.

res <- omicsPCRes_aPC

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$mirna[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$mirna))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from all samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$mirna[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$mirna))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from all samples')


#Nominal p-value plot PC
data_plot <- data.frame(pval=unlist(res$pvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$pvaluesPC[,-c(1,2)]),
                                               each=nrow(res$pvaluesPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values from omicsPC (all samples)')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$qvaluesPC[,-c(1,2)]),
                                               each=nrow(res$qvaluesPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR from omicsPC (all samples)')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)



#Significant genes or miRNA from individual omics (FDR<0.05)

#significant genes
(sig_expr_aPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant miRNA
(sig_mirna_aPC <- rownames(res$mirna)[which(res$mirna[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + miRNA parametric combination (Fisher, FDR<0.05)
sig_expr_mirna_aPC <- res$qvaluesPC[which(res$qvaluesPC[,'Fisher']<0.05),1:2]


#venn diagram
m1 <- list(rna=sig_expr_aPC, rna_fisher=unique(sig_expr_mirna_aPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: PC (all samples)')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(mirna=sig_mirna_aPC, mirna_fisher=unique(sig_expr_mirna_aPC$mirna))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant miRNA identified with individual (mirna) or integration (mirna_fisher)\n omics analysis: PC (all samples)')
grid.newpage()
grid.draw(vennPlot2)



setdiff(sig_expr_mirna_aPC$expr,sig_expr_aPC)[1:50]

setdiff(sig_expr_mirna_aPC$mirna,sig_mirna_aPC)

#Non Parametric Combination
set.seed(12345)

# Setting methods to combine pvalues
combMethods <- c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms <- 1000
# Setting number of cores
numCores <- 4
# Setting omicsNPC to print out the steps that it performs.
verbose <- TRUE

omicsNPCRes_aNPC_rna_mirna <- omicsNPC(dataInput = dataInputExprMirna,
                                       dataMapping = dataMappingExprMirna,
                                       dataTypes = dataTypesExprMirna,
                                       combMethods = combMethods, 
                                       numPerms = numPerms,
                                       numCores = numCores,
                                       verbose = verbose)

res <- omicsNPCRes_aNPC_rna_mirna

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$mirna[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$mirna))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$mirna[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$mirna))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot NPC
data_plot <- data.frame(pval=unlist(res$pvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$pvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$pvaluesNPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsNPC')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$qvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$qvaluesNPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval, fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsNPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)



#Significant genes or miRNA from individual omics (FDR<0.05)

#significant genes
(sig_expr_aNPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant miRNA
(sig_mirna_aNPC <- rownames(res$mirna)[which(res$mirna[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + miRNA non parametric combination (Fisher, FDR<0.05)
sig_expr_mirna_aNPC <- res$qvaluesNPC[which(res$qvaluesNPC[,'Fisher']<0.05),1:2]


#venn diagram
m1 <- list(rna=sig_expr_aNPC, rna_fisher=unique(sig_expr_mirna_aNPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: NPC (all samples)')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(mirna=sig_mirna_aNPC, mirna_fisher=unique(sig_expr_mirna_aNPC$mirna))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant miRNA identified with individual (mirna) or integration (mirna_fisher)\n omics analysis: NPC (all samples)')
grid.newpage()
grid.draw(vennPlot2)



setdiff(sig_expr_mirna_aNPC$expr,sig_expr_aNPC)



setdiff(sig_expr_mirna_aNPC$mirna,sig_mirna_aNPC)

#Integration strategies comparison
#Finally, we can compare the divergences between parametric and non-parametric approaches and between the overlapping or whole cohort escenarios.

#Significant pairs from different approaches

#parmetric combination with overlapping samples (n=397)
oPC <- paste(sig_expr_mirna_cPC$expr,sig_expr_mirna_cPC$mirna, sep=":")
#parmetric combination with all samples (n=466)
aPC <- paste(sig_expr_mirna_aPC$expr,sig_expr_mirna_aPC$mirna, sep=":")
#non parmetric combination with overlapping samples (n=27)
oNPC <- paste(sig_expr_mirna_cNPC$expr,sig_expr_mirna_cNPC$mirna, sep=":")
#non parmetric combination with all samples (n=50)
aNPC <- paste(sig_expr_mirna_aNPC$expr,sig_expr_mirna_aNPC$mirna, sep=":")

res_all <- list(oPC=oPC, aPC=aPC, oNPC=oNPC, aNPC=aNPC)

vennPlot_m1 <- venn.diagram(res_all, NULL, fill=c("green","blue","orange","pink"),
                            cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2,
                            main="Significant pairs for the different escenarios \n FDR (cut-off=0.05)")

grid.newpage()
grid.draw(vennPlot_m1)



#mRNA + Methylation
#Mapping file
#Load the mapping info
#load('methy2GeneMap.RData')
methy2GeneMap <- methy2GeneMap[-40902, ]; #repetition

#restricting mapping info to methylation info
methy2GeneMap_1 <- methy2GeneMap[methy2GeneMap[,1]%in%rownames(logit_met_gbm_corrected),] #from the initial 311,135 methyl features we get 57,645 in common

head(methy2GeneMap_1)

head(expr2GeneMap)

dataMappingExprMet <- combiningMappings(mappings = list(expr = expr2GeneMap, 
                                                        met = methy2GeneMap_1), 
                                        retainAll = TRUE, reference = 'Gene'); 

head(dataMappingExprMet)

#We can see that based on our mRNA and methylation analysis, we have 57,645 pairs in total, it corresponds to 9,620 unique genes (80% from the whole mRNA dataset) and 57,645 unique methylation sites (18% from teh whole methylation dataset).

#Then, we define the relevant clinical variables in each escenario. The variable of interest in the case of GBM is the survivalOtucome. Moreover, Moreover, days_to_birth is included as a confounder.

clinicalVariables <- list(expr = c('days_to_birth', 'survivalOutcome'), 
                          met = c('days_to_birth', 
                                  'survivalOutcome'),
                          mirna = c('days_to_birth', 'survivalOutcome'));
#Overlaping samples
#restrict the datasets to the elements of the data mapping
exprTMP <- rna_gbm_corrected[unique(na.omit(dataMappingExprMet$expr)),
                             ids_rna_met]
#Dimension: 9,620 x 83

metTMP <-logit_met_gbm_corrected[unique(na.omit(dataMappingExprMet$met)),
                                 ids_rna_met]
#Dimension: 57,645 x 83

#specifying the data types.
dataTypesExprMet <- list(ttCoxphContinuous,ttCoxphContinuous)

#preparing the datasets
exprTMP <- createOmicsExpressionSet(Data = as.matrix(exprTMP), 
                                    pData = clinical_gbm[ids_rna_met,
                                                         clinicalVariables[['expr']]])
names(pData(exprTMP)) <- c("age","outcome")

metTMP <- createOmicsExpressionSet(Data = as.matrix(metTMP), 
                                   pData = clinical_gbm[ids_rna_met,
                                                        clinicalVariables[['met']]])
names(pData(metTMP)) <- c("age","outcome")

dataInputExprMet <- list(expr = exprTMP, met = metTMP)
#Parametric Combination
set.seed(12345)
omicsPCRes_cPC <- omicsPC(dataInput = dataInputExprMet,
                          dataMapping = dataMappingExprMet,
                          dataTypes = dataTypesExprMet,
                          verbose = TRUE)

res <- omicsPCRes_cPC

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$met[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$met))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$met[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$met))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot PC
data_plot <- data.frame(pval=unlist(res$pvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$pvaluesPC[,-c(1,2)]),
                                               each=nrow(res$pvaluesPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) + 
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsPC')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$qvaluesPC[,-c(1,2)]),
                                               each=nrow(res$qvaluesPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval,
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


#Now, based on FDR results we can check if there are an increase of statistical power when both omics are jointly analyzed.



#significant genes
(sig_expr_cPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant methyl
(sig_met_cPC <- rownames(res$met)[which(res$met[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + miRNA non parametric combination (Fisher, FDR<0.05)
sig_expr_met_cPC <- res$qvaluesPC[which(res$qvaluesPC[,'Fisher']<0.05),1:2]
#323 pairs

#venn diagram
m1 <- list(rna=sig_expr_cPC, rna_fisher=unique(sig_expr_met_cPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, 
                          main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(met=sig_met_cPC, mirna_fisher=unique(sig_expr_met_cPC$met))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2,
                          main='Significant methylation identified with individual (met) or integration (met_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot2)



setdiff(sig_expr_met_cPC$expr,sig_expr_cPC)[1:50]


setdiff(sig_expr_met_cPC$met,sig_met_cPC)[1:50]

#Non Parametric Combination
set.seed(12345)

# Setting methods to combine pvalues
combMethods <- c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms <- 1000
# Setting number of cores
numCores <- 1
# Setting omicsNPC to print out the steps that it performs.
verbose <- TRUE

omicsNPCRes_cNPC <- omicsNPC(dataInput = dataInputExprMet,
                             dataMapping = dataMappingExprMet,
                             dataTypes = dataTypesExprMet,
                             combMethods = combMethods, 
                             numPerms = numPerms,
                             numCores = numCores,
                             verbose = verbose)

res <- omicsNPCRes_cNPC_rna_met

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$met[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$met))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$met[,"adj.P.Val"]),
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$met))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot NPC
data_plot <- data.frame(pval=unlist(res$pvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$pvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$pvaluesNPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval,
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsNPC')


#FDR plot NPC
data_plot <- data.frame(pval=unlist(res$qvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$qvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$qvaluesNPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval,
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsNPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


Now, based on FDR results we can check if there are an increase of statistical power when both omics are jointly analyzed.

Note: here we consider significant an FDR<0.05.

#Significant genes or miRNA from individual omics (FDR<0.05)

#significant genes
(sig_expr_cNPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])

#significant met
(sig_met_cNPC <- rownames(res$met)[which(res$met[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + methylation non parametric combination (Fisher, FDR<0.05)
sig_expr_met_cNPC <- res$qvaluesNPC[which(res$qvaluesNPC[,'Fisher']<0.05),1:2]


#venn diagram
m1 <- list(rna=sig_expr_cNPC, rna_fisher=unique(sig_expr_met_cNPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, 
                          main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(met=sig_met_cNPC, met_fisher=unique(sig_expr_met_cNPC$met))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2,
                          main='Significant methylation identified with individual (met) or integration (met_fisher)\n omics analysis: PC')
grid.newpage()
grid.draw(vennPlot2)




setdiff(sig_expr_met_cNPC$expr,sig_expr_cNPC)[1:50]

setdiff(sig_expr_met_cNPC$met,sig_met_cNPC)[1:50]

#All samples
#restrict the datasets to the elements of the data mapping
exprTMP <- rna_gbm_corrected[unique(na.omit(dataMappingExprMet$expr)),]; 
#Dimension: 9,620 x 523

metTMP <- logit_met_gbm_corrected[unique(na.omit(dataMappingExprMet$met)),]; 
#Dimension: 57,645 x 95

#specifying the data types.
dataTypesExprMet <- list(ttCoxphContinuous,ttCoxphContinuous)

#preparing the datasets
exprTMP <- createOmicsExpressionSet(Data = as.matrix(exprTMP), 
                                    pData = clinical_gbm[colnames(rna_gbm_corrected),clinicalVariables[['expr']]])

names(pData(exprTMP)) <- c("age","outcome")

metTMP <- createOmicsExpressionSet(Data = as.matrix(metTMP), 
                                   pData = clinical_gbm[colnames(logit_met_gbm_corrected),clinicalVariables[['met']]])

names(pData(metTMP)) <- c("age","outcome")

dataInputExprMet <- list(expr = exprTMP, met = metTMP)
#Parametric Combination
set.seed(12345)
omicsPCRes_aPC <- omicsPC(dataInput = dataInputExprMet,
                          dataMapping = dataMappingExprMet,
                          dataTypes = dataTypesExprMet,
                          verbose = TRUE)

res <- omicsPCRes_aPC

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$met[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$met))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from all samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$met[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$met))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from all samples')


#Nominal p-value plot PC
data_plot <- data.frame(pval=unlist(res$pvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$pvaluesPC[,-c(1,2)]),
                                               each=nrow(res$pvaluesPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values from omicsPC (all samples)')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesPC[,-c(1,2)]), 
                        combination_method=rep(colnames(res$qvaluesPC[,-c(1,2)]),
                                               each=nrow(res$qvaluesPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR from omicsPC (all samples)')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)



#Significant genes or methylation from individual omics (FDR<0.05)

#significant genes
(sig_expr_aPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant methylation sites
(sig_met_aPC <- rownames(res$met)[which(res$met[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + methylation non parametric combination (Fisher, FDR<0.05)
sig_expr_met_aPC <- res$qvaluesPC[which(res$qvaluesPC[,'Fisher']<0.05),1:2]


#venn diagram
m1 <- list(rna=sig_expr_aPC, rna_fisher=unique(sig_expr_met_aPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: PC (all samples)')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(met=sig_met_aPC, met_fisher=unique(sig_expr_met_aPC$met))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant methylation sites identified with individual (met) or integration (met_fisher)\n omics analysis: PC (all samples)')
grid.newpage()
grid.draw(vennPlot2)




setdiff(sig_expr_met_aPC$expr,sig_expr_aPC)

setdiff(sig_expr_met_aPC$met,sig_met_aPC)[1:20]

#Non Parametric Combination
set.seed(12345)

# Setting methods to combine pvalues
combMethods <- c("Fisher", "Liptak", "Tippett")
# Setting number of permutations
numPerms <- 1000
# Setting number of cores
numCores <- 4
# Setting omicsNPC to print out the steps that it performs.
verbose <- TRUE

omicsNPCRes_aNPC_rna_met <- omicsNPC(dataInput = dataInputExprMet,
                                     dataMapping = dataMappingExprMet,
                                     dataTypes = dataTypesExprMet,
                                     combMethods = combMethods, 
                                     numPerms = numPerms,
                                     numCores = numCores,
                                     verbose = verbose)

res <- omicsNPCRes_aNPC_rna_met

#Nominal p-values plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"Pr(>|z|)"],res$met[,"Pr(>|z|)"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)), 
                                rep(names(res)[2], nrow(res$met))))

p1 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics p-value from overlapping samples')

#FDR plot individual omics
data_plot <- data.frame(pval=c(res$expr[,"adj.P.Val"],res$met[,"adj.P.Val"]), 
                        group=c(rep(names(res)[1], nrow(res$expr)),
                                rep(names(res)[2], nrow(res$met))))

p2 <- ggplot(data_plot, aes(x=group, y=pval, fill=group)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of individual omics FDR from overlapping samples')


#Nominal p-value plot NPC
data_plot <- data.frame(pval=unlist(res$pvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$pvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$pvaluesNPC)))

p3 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of p-values in overlapping omicsNPC')


#FDR plot PC
data_plot <- data.frame(pval=unlist(res$qvaluesNPC[,-c(1:4)]), 
                        combination_method=rep(colnames(res$qvaluesNPC[,-c(1:4)])
                                               , each=nrow(res$qvaluesNPC)))

p4 <- ggplot(data_plot, aes(x=combination_method, y=pval, 
                            fill=combination_method)) +
  geom_violin(trim=FALSE)+
  labs(title='Plot of FDR in overlapping omicsNPC')


grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)




#Significant genes or methylation sites from individual omics (FDR<0.05)

#significant genes
(sig_expr_aNPC <- rownames(res$expr)[which(res$expr[,"adj.P.Val"]<0.05)])


#significant methylation sites
(sig_met_aNPC <- rownames(res$met)[which(res$met[,"adj.P.Val"]<0.05)])


#Significant genes from mRNA + met non parametric combination (Fisher, FDR<0.05)
sig_expr_met_aNPC <- res$qvaluesNPC[which(res$qvaluesNPC[,'Fisher']<0.05),1:2]
# 332 pairs

#venn diagram
m1 <- list(rna=sig_expr_aNPC, rna_fisher=unique(sig_expr_met_aNPC$expr))

vennPlot1 <- venn.diagram(m1, NULL, fill=c('red','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant genes identified with individual (rna) or integration (rna_fisher)\n omics analysis: NPC (all samples)')
grid.newpage()
grid.draw(vennPlot1)


m2 <- list(met=sig_met_aNPC, met_fisher=unique(sig_expr_met_aNPC$met))

vennPlot2 <- venn.diagram(m2, NULL, fill=c('green','blue'), 
                          cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, main='Significant methylation sites identified with individual (met) or integration (met_fisher)\n omics analysis: NPC (all samples)')

grid.newpage()
grid.draw(vennPlot2)




setdiff(sig_expr_met_aNPC$expr,sig_expr_aNPC)[1:50]

setdiff(sig_expr_met_aNPC$met,sig_met_aNPC)[1:50]

#Integration strategies comparison

#Significant pairs from different approaches

#parmetric combination with overlapping samples
oPC <- paste(sig_expr_met_cPC$expr,sig_expr_met_cPC$met, sep=":") #323
#parmetric combination with all samples
aPC <- paste(sig_expr_met_aPC$expr,sig_expr_met_aPC$met, sep=":") #207
#non parmetric combination with overlapping samples
oNPC <- paste(sig_expr_met_cNPC$expr,sig_expr_met_cNPC$met, sep=":") #150
#non parmetric combination with all samples
aNPC <- paste(sig_expr_met_aNPC$expr,sig_expr_met_aNPC$met, sep=":") #332

res_all <- list(oPC=oPC, aPC=aPC, oNPC=oNPC, aNPC=aNPC)

vennPlot_m1 <- venn.diagram(res_all, NULL, fill=c("green","blue","orange","pink"), 
                            cat.fontface=4, lty=2, cex=1.2, cat.cex=1.2, 
                            main="Significant pairs for the different escenarios \n FDR (cut-off=0.05)")

grid.newpage()
grid.draw(vennPlot_m1)