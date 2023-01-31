# 
#mirna
# raw_counts <- read.delim("/home/bianca/Downloads/reads_one_per_mirna.csv", sep=",", quote = "")
# precursor<-raw_counts[,4]
# rownames(raw_counts)<-raw_counts[,3]
# exp <- raw_counts[,c(5:10,17:22)]
# colnames(exp)<-str_replace(colnames(exp), "[X][.]", "")
# colnames(exp)<-str_replace(colnames(exp), "[.]", "")
# rownames(exp)<-str_replace(rownames(exp), "[\"]", "")
# rownames(exp)<-str_replace(rownames(exp), "[\"]", "")
# sample_info<-data.frame(colnames(exp))
# sample_info$condition<-c(rep(c(rep(0,3), rep(1,3)),2))
# sample_info$replicate<-as.factor(c(rep(c("1","2","3"), 4)))
# sample_info$times<-c(rep(8,6), rep(24,6))
# sample_info$treatmentL<-c(rep(paste("DMSO","8",sep="_"),3),rep(paste("MKC","8",sep="_"),3),rep(paste("DMSO","24",sep="_"),3),rep(paste("MKC","24",sep="_"),3))
# sample_info$treatment<-as.factor(c(rep(c(rep("DMSO",3), rep("MKC",3)),2)))
# colnames(sample_info)<-c("sample", "condition", "replicate", "times","treatment","drug")
# samples_info_mirna<-sample_info
# write.table(exp,"/home/bianca/mirna_aitor.txt", row.names = T, quote = F,sep="\t")
# mirna<-read.table("/home/bianca/mirna_aitor.txt", sep="\t")
# #normalization
#
# exp1 <- as.data.frame(lapply(mirna, function(x) as.numeric(x)))
# rownames(exp1)<-rownames(mirna) ## && make dge list
# mirna_gbm<- norm_samplesm
# 
# ### rna 
# library("edgeR")
# library("limma")
# library("sva")
# library("biomaRt")
# library(preprocessCore)
# 
# 
# 
# raw_counts <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/countsaitor.txt", sep="", quote = "") 
# raw_counts <- read.delim("/home/bianca/Downloads/countsaitor.txt", sep="", quote = "") 
# rownames(raw_counts)<-raw_counts[,1]
# exp <- raw_counts[,c(7:18)] # first 6 columns is annotation
# sample_info_RNAseq <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/sample_info_RNAseq.txt", sep="") 
# sample_info_RNAseq <- read.delim("/home/bianca/Downloads/sample_info_RNAseq.txt", sep="", row.names = 1) 
# colnames(sample_info_RNAseq)<-c("sample", "replicate", "times","treatment","treat")
# sample_info_RNAseq$condition<-c(rep(0,6), rep(1,6))
# colnames(exp)<-sample_info_RNAseq$sample
# samples_info<-sample_info_RNAseq
# sample_info<-sample_info_RNAseq
# sample_info$short_names<-c("D8R1","D24R3","D24R2","D8R3","D8R2","D24R1","M24R3","M8R2","M8R1","M8R3","M24R2","M24R1")
# colnames(exp)<-sample_info$short_names
# 
# library(biomaRt)
# # #annotation of genes
# httr::set_config(httr::config(ssl_verifypeer = FALSE))
# values=rownames(exp)
# attributes = c("ensembl_gene_id", "hgnc_symbol")
# filters="ensembl_gene_id"
# 
# homo.anno = useEnsembl(biomart = "ensembl", 
#                        dataset = "hsapiens_gene_ensembl", 
#                        mirror = "useast")
# # #listAttributes(homo.anno)
# 
# genes <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
# exp$ens=rownames(exp)
# exp1<-merge(genes, exp, by.x=1,by.y=13)
# library(plyr)
# exp2<-ddply(exp1, .(hgnc_symbol), numcolwise(mean))
# rownames(exp2)=exp2$hgnc_symbol
# exp2<-exp2[,-1]
# 
# y <- DGEList(counts=exp, samples=sample_info_RNAseq , genes=rownames(exp))
# y <- DGEList(counts=exp2, samples=sample_info_RNAseq)
# 
# 
# # REMOVE LOW EXPRESSED GENES
# keep.exprs <- aveLogCPM(y) > 0
# y <- y[keep.exprs, keep.lib.sizes=F]
# 
# # keep.exprs <- filterByExpr(y)
# # y <- y[keep.exprs, keep.lib.sizes=T]
# # 
# # keep.exprs <- rowSums(y$counts==0) < 12
# # y <- y[keep.exprs,, keep.lib.sizes=FALSE]
# 
# 
# #normalization
# 
# function.norm <- function(counts){
#   norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
#   rownames(norm)=rownames(counts)
#   colnames(norm)=colnames(counts)
#   
#   return(norm)
# }
# 
# 
# norm_samples = function.norm(log2(y$counts))
# 
# norm_samples = function.norm(log2(exp2))
# write.table(norm_samples,"/home/bianca/Downloads/rna_normalized_final.txt", sep="\t",quote=F,  row.names=T)
# write.table(sample_info_RNAseq,"/home/bianca/Downloads/rna_normalized_final_sample_info.txt", sep="\t",quote=F,  row.names=T)


######## proteins
# logit_met_gbm<-read.table("/home/bianca/Downloads/protein_matrix_gene_names_zero.txt", sep="\t")
# colnames(logit_met_gbm)<-logit_met_gbm[1,]
# logit_met_gbm<-logit_met_gbm[-1,]
# logit_met_gbm<-logit_met_gbm[!(is.na(logit_met_gbm$`Gene names`) | logit_met_gbm$`Gene names`==""), ]
# logit_met_gbm1<-logit_met_gbm[,c(1:48,64)]
# prot_info<-read.table("/home/bianca/Desktop/proteins_info_last_whole.txt", sep="\t")
# colnames(prot_info)=prot_info[1,]
# prot_info=prot_info[-1,]
# colnames(logit_met_gbm1)=prot_info$short_name
# logit_met_gbm1[, c(1:48)] <- sapply(logit_met_gbm1[, c(1:48)], as.numeric)
# pro<-ddply(logit_met_gbm1, .(Gene_names), numcolwise(mean))
# rownames(pro)=pro$Gene_names
# pro=pro[,-1]
# proteins<-t(apply(pro,1, function(x) tapply(x,colnames(pro),mean)))
# sample_info_proteins<-prot_info[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47),]
# rownames(sample_info_proteins)=sample_info_proteins$short_name 
# sample_info_proteins=sample_info_proteins[colnames(proteins),]
# write.table(proteins,"/home/bianca/Desktop/proteins_final.txt",sep="\t",quote=F)
# write.table(sample_info_proteins,"/home/bianca/Desktop/proteins_sample_info_final.txt",sep="\t",quote=F,row.names = F)


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
mirna_gbm<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final_removed_low_counts.txt",sep="\t")
mirna_gbm_info<-read.table("/home/bianca/Desktop/proteins_inputs/proteins_final_info_correct_order.txt",sep="\t")
proteins=mirna_gbm


rownames(rna_gbm_info)=colnames(rna_gbm)
rna_info_small=as.data.frame(rna_gbm_info[,5])
rownames(rna_info_small)=colnames(rna_gbm)

mm=model.matrix(0 + treat, rna_gbm_info)
B1 <- createOmicsExpressionSet(Data=as.matrix(rna_gbm), pData=rna_info_small, pDataDescr=c("treatment"))
B1 <- createOmicsExpressionSet(Data=as.matrix(rna_gbm), pData=rna_gbm_info)

rna_test=rna_gbm[1:10,]
B1 <- createOmicsExpressionSet(Data=as.matrix(rna_test), pData=rna_info_small, pDataDescr=c("treatment"))

# proteins=mirna_gbm
# y <- DGEList(counts=proteins, samples=sample_info_proteins)
# # REMOVE LOW EXPRESSED GENES
# 
# keep <- rowSums(cpm(y)>3) >= 4
# y <- y[keep, , keep.lib.sizes=FALSE]
# proteins1=y$counts

sample_info_proteins=mirna_gbm_info
# rownames(sample_info_proteins)=sample_info_proteins$V6
# colnames(sample_info_proteins)=sample_info_proteins[1,]
# sample_info_proteins=sample_info_proteins[colnames(proteins),]
# write.table(proteins1,"/home/bianca/Desktop/proteins_final_removed_low_counts.txt",sep="\t",quote=F,row.names=F)
# write.table(sample_info_proteins,"/home/bianca/Desktop/proteins_final_info_correct_order.txt",sep="\t",quote=F,row.names=F)
sample_info_proteins_small=as.data.frame(sample_info_proteins[,5])
rownames(sample_info_proteins_small)=colnames(proteins)

######### time series limma
treatment=sample_info_proteins$joint
levels=c("8HR_DMSO","24HR_DMSO","48HR_DMSO","72HR_DMSO","8HR_MKC","24HR_MKC","48HR_MKC","72HR_MKC")
treatment <- factor(treatment, levels=levels)

## see changes of dmso

#normalize
norm_proteins=function.norm(proteins)
#de time series
design <- model.matrix(~treatment)
fit <- lmFit(norm_proteins,design)
fit <- eBayes(fit,trend=TRUE,robust=FALSE)
a<-capture.output(summary(decideTests(fit[,-1],p.value = 0.1)))
write.table(a,"/home/bianca/Desktop/proteomics_summary_decide_test_time_series.txt",sep="\t",quote=F)
# Intercept = dmso 8
# 24 hours dmso
res24dmso<-topTable(fit,coef=2,n=Inf,p.value = 0.5)
# 48 hours dmso
res48dmso<-topTable(fit,coef=3,n=Inf,p.value = 0.5)
# 72 hours dmso
res72dmso<-topTable(fit,coef=4,n=Inf,p.value = 0.5)
res72dmsoI<-topTable(fit,coef=4,n=Inf,p.value = 0.05)
write.table(res72dmsoI,"/home/bianca/Desktop/proteomics_time_series_72_dmso_whole.txt",sep="\t",quote=F,row.names = T)
write.table(res72dmsoI[,c(1,5)],"/home/bianca/Desktop/proteomics_time_series_72_dmso_bio.txt",sep="\t",quote=F,row.names = T)

# 8 hours mkc
res8mkc<-topTable(fit,coef=5,n=Inf,p.value = 0.5)
# 24 hours mkc
res24mkc<-topTable(fit,coef=6,n=Inf,p.value = 0.5)
# 48 hours mkc
res48mkc<-topTable(fit,coef=7,n=Inf,p.value = 0.5)
# 72 hours mkc
res72mkc<-topTable(fit,coef=8,n=Inf,p.value = 0.5)
# 72 hours mkc important
res72mkcI<-topTable(fit,coef=8,n=Inf,p.value = 0.05)
write.table(res72mkcI,"/home/bianca/Desktop/proteomics_time_series_72_mkc_whole.txt",sep="\t",quote=F,row.names = T)
write.table(res72mkcI[,c(1,5)],"/home/bianca/Desktop/proteomics_time_series_72_mkc_bio.txt",sep="\t",quote=F,row.names = T)

# Block2 - protein expression data
B2 <- createOmicsExpressionSet(Data=as.matrix(proteins), pData=sample_info_proteins_small,pDataDescr=c("treatment"))

samples<-proteins[,c(1:3,10:15,22:24)]
#samples<-proteins1[,c(1:3,10:15,22:24)]
sample_info=sample_info_proteins[c(1:3,10:15,22:24),]
sample_info=as.data.frame(sample_info_proteins_small[c(1:3,10:15,22:24),])
rownames(sample_info)=colnames(samples)
samples_test=samples[1:10,]
B2 <- createOmicsExpressionSet(Data=as.matrix(samples), pData=sample_info)
B2 <- createOmicsExpressionSet(Data=as.matrix(samples_test), pData=sample_info,pDataDescr=c("treatment"))

ms <- modelSelection(Input=list(B1, B2), Rmax=6, fac.sel="single%",
                     varthreshold=0.03, center=TRUE, scale=TRUE, 
                     weight=TRUE, plot_common=FALSE, plot_dist=FALSE)

ms <- modelSelection(Input=list(B1, B2), Rmax=1000, fac.sel="single%",
                     varthreshold=0.03, center=TRUE, scale=TRUE, 
                     weight=TRUE, plot_common=TRUE, plot_dist=TRUE)


discoRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("rna", "proteins"),
                              method="DISCOSCA", Rcommon=3, Rspecific=c(4, 8),
                              center=TRUE, scale=TRUE, weight=TRUE)
jiveRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("rna", "proteins"),
                             method="JIVE", Rcommon=3, Rspecific=c(4, 8),
                             center=TRUE, scale=TRUE, weight=TRUE)
o2plsRes <- omicsCompAnalysis(Input=list(B1, B2),Names=c("rna", "proteins"),
                              method="O2PLS", Rcommon=3, Rspecific=c(4, 8),
                              center=TRUE, scale=TRUE, weight=TRUE)

dataTypes <- c("count", "count")

combMethods <- c("Fisher", "Liptak", "Tippett")
combMethods <- c("Tippett")

numPerms <- 1000

numCores <- 1

verbose <- TRUE


results <- omicsNPC(dataInput=list(B1, B2), dataTypes=dataTypes,
                    combMethods=combMethods, numPerms=numPerms,
                    numCores=numCores, verbose=verbose)

