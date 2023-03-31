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

option_list = list(
  make_option(
    "--samples",
    action = "store",
    default = "/home/ben/Desktop/pacli/tapsi_star/t8_star_tab/PD8.txt",
    type = 'character',
    help = 'Path to samples files, for multiple files add *.txt'
  ),
  make_option(
    "--samplesInfo",
    action = "store",
    default = "/home/ben/Desktop/pacli/tapsi_star/samples_info_8.csv",
    type = 'character',
    help = 'Path to samplesInfo files, for multiple files add *.txt'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))
print("Importing")

allsamples <- lapply(Sys.glob(opt$samples), read.table)
allsamples_info <- lapply(Sys.glob(opt$samplesInfo), read.table)
allnames<-lapply(Sys.glob(opt$samples), ohmyname)

#########################################################
# keep name as variable
ohmyname<-function(mrna){
  whoa<-str_extract(mrna, regex("[^/]+$*[\\.]"))
  whoa2<-sub("\\.","",whoa)
  return(whoa2)
}

#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  return(norm)
}


##3 paclitaxel
counts1<-read.table("/home/ben/Desktop/pacli/tapsi_star/t8_star_tab/PD8.txt")
colnames(counts1)=counts1[1,]
counts1=counts1[-1,]
counts2<-read.table("/home/ben/Desktop/pacli/tapsi_star/t24_star_tab/PD24.txt", sep=",")
colnames(counts2)=counts2[1,]
counts2=counts2[-1,]

counts<-cbind(counts1,counts2)
rownames(counts)=counts[,1]
counts=counts[,-c(1,8)]

write.table(counts, "/home/ben/Desktop/pacli/tapsi_star/all_together_pac.txt", sep="\t", quote=F, row.names = T)

si8<-read.table("/home/ben/Desktop/pacli/tapsi_star/samples_info_8.csv")
si24<-read.table("/home/ben/Desktop/pacli/tapsi_star/samples_info_24.csv")

si<-rbind(si8,si24)
colnames(si)=si[1,]
si=si[-c(1,8),]
si$groups<-paste(si$treatment, si$hours, sep="")
write.table(si, "/home/ben/Desktop/pacli/tapsi_star/all_together_pac_samplesInfo.txt", sep="\t", quote=F, row.names = T)

counts<-read.table("/home/ben/Desktop/pacli/all_together_pac.txt")
si<-read.table("/home/ben/Desktop/pacli/all_together_pac_samplesInfo.txt")

#################################################################################################################


# create DGEobject 
# for( i in allsamples){
#   for( x in allsamples_info){
#     y <- list(DGEList(counts=as.data.frame(allsamples[i]), samples=allsamples_info[x], genes=rownames(as.data.frame(allsamples[i]))))
#   }
# }


counts <- mutate_all(counts, function(x) as.numeric(as.character(x)))
y <- DGEList(counts=counts, samples=si, genes=rownames(counts))

names(y)

y$counts<-cpm(y, log=TRUE)
y$counts<-cpm(y)
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

# REMOVE LOW EXPRESSED GENES
table(rowSums(y$counts==0)==12)

keep.exprs <- filterByExpr(y,design=design, min.count=3)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]


L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)


####### normalization
norm_samples<-function.norm(y$counts)

#### or

y <- calcNormFactors(y, method = "upperquartile")
y$samples$norm.factors



lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")



#################### visualization
x2 <- y
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm2 <- cpm(x2, log=TRUE)
boxplot(lcpm2, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#############################


############################ MDS Unsupervised clustering of samples
lcpm3 <- cpm(y, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm3, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm3, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

######## or
glMDSPlot(lcpm3, labels=paste(group, lane, sep="_"), 
          groups=y$samples[,c(2,5)], launch=FALSE)


####################################################################


norm_samples<-function.norm(y$counts)
#### OR



# batch effect corrections

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



################
#batch effect correction
colnames(norm_samples)=si$sample
comsvaresres<-comsva(norm_samples, si , treatment, replicate)
PreparePlotpca(comsvaresres[[1]], si, si$treatment, si$replicate, cbp_12)
svacomrere<-svacom(norm_samples, si , treatment, replicate)
PreparePlotpca(svacomrere[[1]], si, si$treatment, si$replicate, cbp_12)
svaresres<-svAll(norm_samples, si , treatment, replicate)
#PreparePlotpca(svaresres[[1]], si, si$treatment, si$replicate, cbp_12)
comresres<-com(norm_samples, si , replicate)
PreparePlotpca(comresres, si, si$treatment, si$replicate, cbp_12)


## batch effect pacli
sv<-list(comsva(y$counts, si, treatment, replicate))

comsvaresres<-comsva(y$counts, si , treatment, replicate)


############################ DIFFERENTIAL EXPRESSION ANALYSIS EdgeR


## mrna

design <- model.matrix(~0 + treat, si)

## pacli
design <- model.matrix(~0 + groups, si)

## design with batch effect
## pacli
design1 <- cbind(design,comsvaresres[[2]][[1]])

### mrna
design<-cbind(design, comsvaresres[[2]][[1]])

colnames(design)<-make.names(colnames(design))
colnames(design1)<-make.names(colnames(design1))

y1<-DGEList(counts=comsvaresres[[1]], samples=si ,genes=row.names(comsvaresres[[1]]) )

## mrna
y<-DGEList(counts=norm_samples, samples=si ,genes=row.names(norm_samples) )
y <- estimateDisp(y, design)
y <- estimateDisp(y, design, robust = T)
y <- estimateGLMTagwiseDisp(y, design)

##
y2 <- estimateDisp(y, design1, robust = T)
y2 <- estimateGLMCommonDisp(y,design1)
y2 <- estimateGLMTagwiseDisp(y, design1)

##
y3 <- estimateDisp(y1, design, robust = T)
y3 <- estimateGLMCommonDisp(y1,design)
y3 <- estimateGLMTagwiseDisp(y1, design)


#fit <- glmFit(y, design)
fit <- glmQLFit(y, design)

## 
fit2 <- glmQLFit(y2, design1)


## mrna
contr.matrix <- makeContrasts(
  mkcdmso8 = treatMKC_8-treatDMSO_8, 
  mkcdmso24 = treatMKC_24 - treatDMSO_24,
  levels = colnames(design))

contr.matrix <- makeContrasts(
  mkcdmso8 = treatMKC_8-treatDMSO_8,
  levels = colnames(design))


## pacli
contr.matrix <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  pacdmso24 = groupspac24 - groupsdmso24,
  levels = colnames(design))

contr.matrix1 <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  pacdmso24 = groupspac24 - groupsdmso24,
  levels = colnames(design1))

contr.matrix8 <- makeContrasts(
  pacdmso8 = groupspac8-groupsdmso8, 
  levels = colnames(design))

contr.matrix24 <- makeContrasts(
  pacdmso24 = groupspac24-groupsdmso24, 
  levels = colnames(design))


#####################3
## pacli
pacdmsores <- glmQLFTest(fit, contrast=contr.matrix)
pacdmsores2 <- glmQLFTest(fit2, contrast=contr.matrix1)
fres<-topTags(pacdmsores, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
fres5<-topTags(pacdmsores, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

fres2<-topTags(pacdmsores2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1)
saveRDS(pacdmsores,"/home/ben/Desktop/pacli/tapsi_star/pacdmsores")
fres<-topTags(pacdmsores, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
fres
write.table(fres, "/home/ben/Desktop/pacli/limma_res_all.txt", sep="\t", quote=F, row.names = F)


f8<-fres[,c(2,7)]
f24<-fres[,c(3,7)]
write.table(f8, "/home/ben/Desktop/pacli/limma_res_8bio.txt", sep="\t", quote=F, row.names = T)
write.table(f24, "/home/ben/Desktop/pacli/limma_res_24bio.txt", sep="\t", quote=F, row.names = T)
write.table(f8, "/home/ben/Desktop/pacli/limma_res_8integrate.txt", sep="\t", quote=F, row.names = T)
write.table(f24, "/home/ben/Desktop/pacli/limma_res_24_integrate.txt", sep="\t", quote=F, row.names = T)

pacdmsores8 <- glmQLFTest(fit, contrast=contr.matrix8)
saveRDS(pacdmsores8,"/home/ben/Desktop/pacli/tapsi_star/pacdmsores8")
fres8<-topTags(pacdmsores8, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
fres8

pacdmsores24 <- glmQLFTest(fit, contrast=contr.matrix24)
saveRDS(pacdmsores24,"/home/ben/Desktop/pacli/tapsi_star/pacdmsores24")
fres24<-topTags(pacdmsores24, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
fres24

## mrna
mkcdmsores <- glmQLFTest(fit, contrast=contr.matrix)
fres<-topTags(mkcdmsores, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
fresFdr<-topTags(mkcdmsores, n = Inf, adjust.method = "fdr", sort.by = "PValue", p.value = 0.05)


#pacdmso8res <- glmQLFTest(fit, coef=2)
#topTags(pacdmso8res)

### orrr likehood test
fit <- glmFit(y, design)

pacdmso8res <- glmLRT(fit, coef=2)
topTags(pacdmso8res)



par(mfrow=c(1,2))
v <- voom(y, design, plot=TRUE)
v2 <- voom(y, design1, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit2 <- lmFit(v2, design1)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
vfit2 <- contrasts.fit(vfit2, contrasts=contr.matrix1)

efit <- eBayes(vfit)
efit2 <- eBayes(vfit2)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))


tfit <- treat(vfit, lfc=0.5)
tfit2 <- treat(vfit2, lfc=0.5)
tfit2 <- treat(vfit2)
tfit <- treat(vfit)
dt <- decideTests(tfit, p.value = 1)
dt2 <- decideTests(tfit2, p.value = 1)
dt2 <- decideTests(tfit2, p.value = 0.05)
summary(dt)
summary(dt2)


###########################################################