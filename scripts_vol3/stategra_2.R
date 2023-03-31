library(STATegRa)
library(parallel)

# oct 22
mirna_gbm_norm<-readRDS("/home/ben/Downloads/de_p_8_normal_values")
rna_gbm_norm<-readRDS("/home/ben/Downloads/de_r_8_normal_values")
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
rownames(samples_info)=colnames(mirna_gbm_norm)
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
colnames(samples_info)=c("condition")
rownames(samples_info)=colnames(mirna_gbm_norm)
mda<-list(mirna=mirna_gbm_norm,rna=rna_gbm_norm)

B1 <- createOmicsExpressionSet(Data=rna_gbm_norm, pData=samples_info,
                               pDataDescr=c("classname"))

B2 <- createOmicsExpressionSet(Data=as.matrix(mirna_gbm_norm), pData=samples_info,
                               pDataDescr=c("classname"))


# # Block1 - gene expression data
# B1 <- createOmicsExpressionSet(Data=rna_gbm_norm, pData=samples_info,
#                                pDataDescr=c("classname"))
# B1 <- createOmicsExpressionSet(Data=rna_gbm_norm, pData=samples_info)
# 
# # sep 22
# 
# rna<-readRDS("/home/bianca/Desktop/rna_filt_norm.RDS")
# r1=rna[,c(1,2,3,7,8,9)]
# rn<-as.matrix(log2(r1+1))
# samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
# rownames(samples_info)=colnames(rn)
# samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
# colnames(samples_info)=c("condition")
# rownames(samples_info)=colnames(r1)
# 
# B1 <- createOmicsExpressionSet(Data=rn, pData=samples_info,pDataDescr=c("condition"))
# 
# ###
# rownames(sample_info_proteins)=sample_info_proteins$common_name
# 
# proteins<-t(apply(met_gbm_norm,1, function(x) tapply(x,colnames(met_gbm_norm),mean)))
# sample_info_proteins1<-sample_info_proteins[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47),]
# rownames(sample_info_proteins1)=sample_info_proteins1$common_name
# sample_info_proteins1=sample_info_proteins1[colnames(proteins),]
# # Block2 - protein expression data
# 
# ## sep 22
# prot<-readRDS("/home/bianca/Desktop/prot.RDS")
# x_b=rownames(prot)
# prot1=as.data.frame(sapply(prot, function(x) as.numeric(x)))
# rownames(prot)=x_b
# pr1=prot[,c(1,2,3,7,8,9)]
# pr1<-as.matrix(log2(pr1+1))
# B2 <- createOmicsExpressionSet(Data=pr1, pData=samples_info,pDataDescr=c("condition"))
# 
# B2 <- createOmicsExpressionSet(Data=Block1.PCA, pData=ed.PCA)
# B2 <- createOmicsExpressionSet(Data=Block2.PCA, pData=ed.PCA)

ms <- modelSelection(Input=list(B1, B2), Rmax=3, fac.sel="single%",
                     varthreshold=0.03, center=TRUE, scale=TRUE, 
                     weight=TRUE, plot_common=FALSE, plot_dist=FALSE)
# ms <- modelSelection(Input=list(B1, B2),center=TRUE, scale=TRUE, 
#                      weight=TRUE, plot_common=FALSE, plot_dist=FALSE)

grid.arrange(ms$common$pssq, ms$common$pratios, ncol=2)

discoRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "proteins"),
                              method="DISCOSCA", Rcommon=3, Rspecific=c(2, 2),
                              center=F, scale=T, weight=TRUE)
jiveRes <- omicsCompAnalysis(Input=list(B1, B2), Names=c("expr", "proteins"),
                             method="JIVE", Rcommon=3, Rspecific=c(2, 2),
                             center=TRUE, scale=TRUE, weight=TRUE)
o2plsRes <- omicsCompAnalysis(Input=list(B1, B2),Names=c("expr", "proteins"),
                              method="O2PLS", Rcommon=3, Rspecific=c(2, 2),
                              center=TRUE, scale=TRUE, weight=TRUE)

dataTypes=c("count", "count")
dataTypes=c("count", "count", "continuous")
combMethods=c("Fisher", "Liptak","Tippett")
mda<-list(rna=B1, prot=B2)
#mda<-TCGA_BRCA_Data

results<-omicsNPC(dataInput = mda, dataTypes=dataTypes, combMethods = combMethods, numPerms=10, numCores=4, verbose=TRUE)
results<-omicsPC(dataInput = mda, dataTypes=dataTypes, verbose=TRUE) # ******** DOESN'T EXIST!!!!!!!!!!!!
