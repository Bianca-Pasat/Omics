library(STATegRa)
library(MASS)
library(gridExtra)
library(devtools)
library('RpESCA')
library(preprocessCore)
library(plyr)
library(RegularizedSCA)


rna_gbm<-read.table("/home/bianca/Downloads/rna_normalized_final.txt",sep="\t")
rna_gbm_info<-read.table("/home/bianca/Downloads/rna_normalized_final_sample_info.txt",sep="\t")
mirna_gbm<-read.table("/home/bianca/Desktop/proteins_final.txt",sep="\t")
mirna_gbm_info<-read.table("/home/bianca/Desktop/proteins_sample_info_final.txt",sep="\t")

rna_gbm<-rn
mirna_gbm<-pr1

#Calculate Forbenious normalization
frobenius_rna <- norm(as.matrix(rna_gbm), type="F")
frobenius_mirna <- norm(as.matrix(mirna_gbm), type="F")
frobenius_met <- norm(as.matrix(logit_met_gbm), type="F")

#Mean-centring and division by Frobenius normalization factor
rna_gbm_norm <- t(scale(t(rna_gbm), scale=FALSE))/frobenius_rna
mirna_gbm_norm <- t(scale(t(mirna_gbm), scale=FALSE))/frobenius_mirna
met_gbm_norm <- t(scale(t(logit_met_gbm), scale=FALSE))/frobenius_met

##Mean-centring without Frobenius normalization
rna_gbm_norm2 <- t(scale(t(rna_gbm), scale=FALSE))
mirna_gbm_norm2 <- t(scale(t(mirna_gbm), scale=FALSE))
met_gbm_norm2 <- t(scale(t(logit_met_gbm), scale=FALSE))


#Common samples between mRNA and mirna
ids_rna_mirna <- intersect(colnames(rna_gbm), colnames(mirna_gbm)) 
length(ids_rna_mirna)

#Common samples between mRNA and proteomics
ids_rna_met <- intersect(colnames(rna_gbm), colnames(logit_met_gbm))
length(ids_rna_met)
                   
library(r.jive)

#Using r.jive

####################################################
## mRNA + miRNA (already centered and normalized) ##
####################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_mirna],
             miRNA=mirna_gbm_norm[,ids_rna_mirna])
rjive_rna_mirna <- jive(Data, scale = FALSE, center = FALSE) 


##########################################################
## mRNA + methylation (already centered and normalized) ##
##########################################################
Data <- list(mRNA=rna_gbm_norm[,ids_rna_met],
             met=met_gbm_norm[,ids_rna_met])
rjive_rna_met <- jive(Data, scale = FALSE, center = FALSE)
# Load require packages
library(RegularizedSCA)

#Using pca-gca

####################################################
## mRNA + miRNA (already centered and normalized) ##
####################################################
data <- cbind(t(rna_gbm_norm[,ids_rna_mirna]), 
              t(mirna_gbm_norm[,ids_rna_mirna]))

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
  alpha_est <- alpha_estimation(X,K=3,opts=opts)
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
sel1 <- which(pESCA_L2$varExpPCs[1,]>5) #3
sel2 <- which(pESCA_L2$varExpPCs[2,]>5) #2
intersect(sel1,sel2)
# 1 common components
# 2 dist mRNA + 1 miRNA
# Total: 4 components


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


#Cut-off 5%
sel1 <- which(pESCA_L2$varExpPCs[1,]>5) #3
sel2 <- which(pESCA_L2$varExpPCs[2,]>5) #3
intersect(sel1,sel2)


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

p <- ggplot(df,aes_string("PC2","PC7",color="batch"))
p + geom_point(size=5, shape=20) + ggtitle("PCA by batch") + theme(legend.position="bottom")

df <- data.frame(common_pc,dist_rna_pc,dist_met_pc)
colnames(df) <- c(paste("common_",colnames(common_pc),sep=""),
                  paste("dist_rna_",colnames(dist_rna_pc),sep=""),
                  paste("dist_met_",colnames(dist_met_pc),sep=""))

pcaCorrelation(clinical=clin, clin_var, df, maxlogPvalue =10, maxPC=ncol(df), pc.names=TRUE)

# Load require package
library(sva)
#Before to apply ComBat, it’s important to note that there are two batches with only one sample

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

pcaCorrelation(clinical=clinical_gbm[colnames(logit_met_gbm_corrected),],
               clin_var=clin_var, pca_data=pca_met$x, maxlogPvalue = 10)

