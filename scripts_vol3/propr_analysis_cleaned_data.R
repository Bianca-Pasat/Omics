load("/home/ben/Desktop/hell.RData")



X_info<-read.table("/home/ben/Downloads/proteomics_all/proteins_inputs/proteins_final_info_correct_order.txt",sep="\t")
colnames(X_info)=X_info[1,]
X_info=X_info[-1,]
X_info$sample=colnames(X)
colnames(X_info)[2]<-("times")
colnames(X_info)[5]<-("treatment")
#colnames(X_info)[6]<-("sample")
X_infoB<-X_info

X<-read.table("/home/ben/Downloads/proteomics_all/proteins_inputs/proteins_final.txt",sep="\t")
x_b=X
X=as.data.frame(lapply(X, function(x) as.numeric(x)))
rownames(X)=rownames(x_b)


## exclude low counts
#keep=apply(X,2,function(x) sum(x>=2) >= 3)
keep=X[!(rowSums(X==0) >=2),]
X1=X[rownames(keep),]
saveRDS(X1,"/home/ben/Desktop/proteomics_cleaned_bianca/removed_low_counts_proteins")
X1<-function.norm(as.matrix(log2(X1)))
saveRDS(X1,"/home/ben/Desktop/proteomics_cleaned_bianca/removed_low_counts_proteins_QUANTILE_NORM_log2")
################################## PCA AND BATCH EFFECT
### 8 hours 
X1=X1[,c(10:12,22:24)]
X_info=X_infoB[c(10:12,22:24),]

PreparePlotpca(X1, X_info, X_info$treatment, X_info$replicate, cbp_12 )

comsvaRe<-comsva(X1, X_info, treatment, replicate)
PreparePlotpca(comsvaRe, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)
PreparePlotpca(comsvaRe, X_info,X_info$condition,as.factor(X_info$times),values=cbp_12)

svacomre<-svacom(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(svacomre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

comre<-com(X1, X_info, replicate)
PreparePlotpca(comre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

savre<-svAll(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(savre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

clean_prot8<-comsvaRe

# 24 hours
X1=X[rownames(keep),]
X1<-function.norm(as.matrix(log2(X1+1)))
X1=X1[,c(1:3,13:15)]
X_info=X_infoB[c(1:3,13:15),]
PreparePlotpca(X1, X_info, X_info$treatment, X_info$replicate, cbp_12 )

comsvaRe<-comsva(X1, X_info, treatment, replicate)
PreparePlotpca(comsvaRe, X_info,X_info$condition,as.factor(X_info$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(svacomre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

comre<-com(X1, X_info, replicate)
PreparePlotpca(comre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

savre<-svAll(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(savre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

clean_prot24<-comsvaRe

# 48 hours
X1=X[rownames(keep),]
X1<-function.norm(as.matrix(log2(X1+1)))
X1=X1[,c(4:6,16:18)]
X_info=X_infoB[c(4:6,16:18),]
PreparePlotpca(X1, X_info, X_info$treatment, X_info$replicate, cbp_12 )

comsvaRe<-comsva(X1, X_info, treatment, replicate)
PreparePlotpca(comsvaRe, X_info,X_info$condition,as.factor(X_info$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(svacomre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

comre<-com(X1, X_info, replicate)
PreparePlotpca(comre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

savre<-svAll(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(savre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

clean_prot48<-comsvaRe
clean_prot48_b<-svacomre

# 72 hours
X1=X[rownames(keep),]
X1<-function.norm(as.matrix(log2(X1+1)))

X1=X1[,c(7,8,9,19,20,21)]
X_info=X_infoB[c(7:9,19:21),]
PreparePlotpca(X1, X_info, X_info$treatment, X_info$replicate, cbp_12 )

comsvaRe<-comsva(X1, X_info, treatment, replicate)
PreparePlotpca(comsvaRe, X_info,X_info$condition,as.factor(X_info$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(svacomre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

comre<-com(X1, X_info, replicate)
PreparePlotpca(comre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

savre<-svAll(as.matrix(X1), X_info, treatment,replicate)
PreparePlotpca(savre, X_info,X_info$treatment,as.factor(X_info$replicate),values=cbp_12)

clean_prot72<-comsvaRe

saveRDS(clean_prot8,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot8")
saveRDS(clean_prot24,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot24")
saveRDS(clean_prot48,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot48")
saveRDS(clean_prot48_b,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot48_B")
saveRDS(clean_prot72,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot72")

#### RANKPROD LOAD THE RDS AND RERUN RANKPROD FOR EACH TIME POINT - DONE
is_logged=TRUE
data_size=FALSE
is_normalized=TRUE
cl=c(rep(1,3),rep(0,3))

clean_prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot8")
clean_prot24<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot24")
clean_prot48<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot48")
clean_prot72<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot72")

# 8
samples=clean_prot8

out_dir="/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8.txt"
upsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8_DEregulatedBio.txt"

# 24
samples=clean_prot24
out_dir="/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_24.txt"
upsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_24_upregulatedBio.txt"
downsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_24_downregulatedBio.txt"
allsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_24_DEregulatedBio.txt"

# 48
samples=clean_prot48
out_dir="/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_48.txt"
upsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_48_upregulatedBio.txt"
downsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_48_downregulatedBio.txt"
allsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_48_DEregulatedBio.txt"

# 72
samples=clean_prot72
out_dir="/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_72.txt"
upsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_72_upregulatedBio.txt"
downsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_72_downregulatedBio.txt"
allsave= "/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_72_DEregulatedBio.txt"


### run analysis
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

########## keep DE proteins with clean matrix
de8<-read.table("/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
de8<-rownames(de8)
de24<-read.table("/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_24_DEregulatedBio.txt",sep="\t")
de24<-rownames(de24)
de48<-read.table("/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_48_DEregulatedBio.txt",sep="\t")
de48<-rownames(de48)
de72<-read.table("/home/ben/Desktop/proteomics_cleaned_bianca/rp_clean_prot_mkc_vs_dmso_72_DEregulatedBio.txt",sep="\t")
de72<-rownames(de72)

cleaned_de_prot8<-clean_prot8[rownames(clean_prot8) %in% de8, ]
cleaned_de_prot24<-clean_prot24[rownames(clean_prot24) %in% de24, ]
cleaned_de_prot48<-clean_prot48[rownames(clean_prot48) %in% de48, ]
cleaned_de_prot72<-clean_prot72[rownames(clean_prot72) %in% de72, ]

saveRDS(cleaned_de_prot8,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot8")
saveRDS(cleaned_de_prot24,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot24")
saveRDS(cleaned_de_prot48,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot48")
saveRDS(cleaned_de_prot72,"/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot72")

############### keep filtered and normalized log2 matrix with de #######################
fitnorm<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/removed_low_counts_proteins_QUANTILE_NORM_log2")

fitnorm8<-fitnorm[rownames(fitnorm) %in% de8, ]
fitnorm8=fitnorm8[,c(10:12,22:24)]
fitnorm24<-fitnorm[rownames(fitnorm) %in% de24, ]
fitnorm24=fitnorm24[,c(1:3,13:15)]
fitnorm48<-fitnorm[rownames(fitnorm) %in% de48, ]
fitnorm48=fitnorm48[,c(4:6,16:18)]
fitnorm72<-fitnorm[rownames(fitnorm) %in% de72, ]
fitnorm72=fitnorm72[,c(7,8,9,19,20,21)]

saveRDS(fitnorm8,"/home/ben/Desktop/proteomics_cleaned_bianca/fitnorm_de_prot8")
saveRDS(fitnorm24,"/home/ben/Desktop/proteomics_cleaned_bianca/fitnorm_de_prot24")
saveRDS(fitnorm48,"/home/ben/Desktop/proteomics_cleaned_bianca/fitnorm_de_prot48")
saveRDS(fitnorm72,"/home/ben/Desktop/proteomics_cleaned_bianca/fitnorm_de_prot72")


######## propr analysis

X=t(X)
X1=X[,keep]
#8 hours
X1=X1[c(10:12,22:24),]
X_info=X_info[c(10:12,22:24),]
# 24 hours
X1=X1[c(1:3,13:15),]
X_info=X_info[c(1:3,13:15),]
# 48 hours
X1=X1[c(4:6,16:18),]
X_info=X_info[c(4:6,16:18),]
# 72 hours
X1=X1[c(7:9,19:21),]
X_info=X_info[c(7:9,19:21),]

# #### PCA RAW
# X2=t(X1)
# 
# PreparePlotpca2(X2, X_info, joint, replicate, cbp_12  )

rho2= propr(X1, metric="rho")
# choose cut off for FDR to be below 0.05
updateCutoffs(rho2, cutoff = seq(.989,.990,.001),ncores=6)
updateCutoffs(rho2, cutoff = seq(.995,.999,.001),ncores=6)

# subset with that cutoff
best=rho2[">",.990]
best=rho2[">",.999] # 24 hours
best=rho2[">",.995] # 72 hours

## 
plot(best)
dendrogram(best)
snapshot(best)

### plots to identify which cluster you should use for pca (Choose the )
best=simplify(best)

k=2
pdf("/home/bianca/48clust")
clusts=prism(best, k=k)
clusts=prism(best, k=k, prompt=F)
dev.off()
clusts=prism(best,k=5)
clusts=bokeh(best, k=k)

## you choose the best cluster

cll=1
sub=subset(best, select= (clusts==cll))
pdf("/home/bianca/48pcacl1")
pca(sub, group= X_info$condition)
dev.off()
transcripts=colnames(sub@logratio)
write.table(transcripts,"/home/bianca/Desktop/transcripts_8.txt", sep="\t",quote=F, row.names = TRUE)
write.table(transcripts,"/home/bianca/Desktop/transcripts_24.txt", sep="\t",quote=F, row.names = TRUE)
write.table(transcripts,"/home/bianca/Desktop/transcripts_48.txt", sep="\t",quote=F, row.names = TRUE)
write.table(transcripts,"/home/bianca/Desktop/transcripts_72.txt", sep="\t",quote=F, row.names = TRUE)

##
sub72=sub
cytescape(best[">", .998])
#cytescape(best.1[">", .95], prompt = FALSE, col1 = TTpos, col2 = TTneg, d3 = TRUE)
##########################################################


#################      Differential proportionality ##########################

## actual analysis
pd=propd(X1, X_info$joint, alpha=0.01,weighted = TRUE, p =100)
theta_d=setDisjointed(pd)
theta_e=setEmergent(pd)
pd.nn=updateF(theta_d, moderated = TRUE, ivar="clr")
tab<-getResults(pd.nn)
t1<-tab[tab$FDR < 0.05,]
pairs=c(t1$Partner, t1$Pair)

# compare with rankprod
de8=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
de24=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_24_DEregulatedBio.txt",sep="\t")
de48=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_48_DEregulatedBio.txt",sep="\t")
de72=read.table("/home/bianca/Desktop/proteomics_de_limma_normal_rankprod_results/rp_prot_mkc_vs_dmso_72_DEregulatedBio.txt",sep="\t")

de8$time=rep("8",length(rownames(de8)))
de24$time=rep("24",length(rownames(de24)))
de48$time=rep("48",length(rownames(de48)))
de72$time=rep("72",length(rownames(de72)))

de_all=rbind(de8,de24,de48,de72)
isoforms=c("SLMAP","REPIN1","HLA-C","HLAC","OSBPL1A","EGLN1","IMMT","NAA16","SUOX","EIF4A3","TUBB6")
d_is_ex_usage=isoforms[isoforms %in% rownames(de_all)]

import=unique(pairs[pairs %in% trans$V1])
import_rank=unique(pairs[pairs %in% rownames(de_prot_8)])
import_rank=unique(pairs[pairs %in% rownames(de8)])
import_rank=unique(pairs[pairs %in% rownames(de24)])
import_rank=unique(pairs[pairs %in% rownames(de48)])
import_rank=unique(pairs[pairs %in% rownames(de72)])

import_rank2=de8[rownames(de8) %in% pairs,]
import_rank2=de24[rownames(de24) %in% pairs,]
import_rank2=de48[rownames(de48) %in% pairs ,]
import_rank2=de72[rownames(de72) %in% pairs ,]


import_trans=trans[rownames(de_prot_8) %in% trans]
super_import=unique(pairs[pairs %in% import_trans])

write.table(import_trans,"/home/bianca/Desktop/common_rank_propr_8.txt", sep="\t", quote=F, row.names = F)
write.table(super_import,"/home/bianca/Desktop/important_and_common_rank_propr_8.txt", sep="\t", quote=F, row.names = F)

write.table(import_rank,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_8.txt", sep="\t", quote=F, row.names = F)
write.table(import_rank,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_24.txt", sep="\t", quote=F, row.names = F)
write.table(import_rank,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_48.txt", sep="\t", quote=F, row.names = F)
write.table(import_rank,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_72.txt", sep="\t", quote=F, row.names = F)

write.table(import_rank2,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_8_FC.txt", sep="\t", quote=F, row.names = T)
write.table(import_rank2,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_24_FC.txt", sep="\t", quote=F, row.names = T)
write.table(import_rank2,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_48_FC.txt", sep="\t", quote=F, row.names = T)
write.table(import_rank2,"/home/bianca/Desktop/deModules_with_propr_and_rankprod_72_FC.txt", sep="\t", quote=F, row.names = T)


############################################################## notes notes notes notes ################################################################3
####################################################### NOTES
####################################### other metrics of proportionality
phi=propr(X, metric="phi", symmetrize = TRUE) # *** see if you change that what happens(syymetrize)
rho=propr(X, metric="rho", ivar=0) ## same as before but on whole
phs=propr(X, metric="phs", ivar=0)

## evaluating proportionality
updateCutoffs(phi, cutoff = seq(.05,.95,.3), ncores=5)
updateCutoffs(rho, cutoff = seq(.05,.99,.3),ncores=6)
updateCutoffs(phs, cutoff = seq(.05,.95,.3), ncores=7)

# rho2= propr(X1, metric="rho", select=keep)
# phi2= propr(X, metric="phi", select=keep) ## here you can put on select lets say the DE proteins.
# best=phi2[">",.9]
# 
# best=rho2[">",.9]
#disjoint proportionality : slope changes
tab=getResults(theta_d)
plot(theta_d@counts[,7], theta_d@counts[,42], col=ifelse(theta_d@group=="DMSO","red","blue"))
grp1=theta_d@group=="DMSO"
grp2=theta_d@group=="MKC"
abline(a=0,b=theta_d@counts[grp1,7]/theta_d@counts[grp1,42], col="red")
abline(a=0,b=theta_d@counts[grp2,7]/theta_d@counts[grp2,42], col="red")
plot(theta_d@counts[,42]/theta_d@counts[,7], col=ifelse(theta_d@group=="DMSO","red","blue"))

# emergent proportionality : strength changes
tab=getResults(theta_e)
plot(theta_e@counts[,7], theta_e@counts[,42], col=ifelse(theta_e@group=="DMSO","red","blue"))
grp1=theta_e@group=="DMSO"
grp2=theta_e@group=="MKC"
abline(a=0,b=theta_e@counts[grp1,7]/theta_e@counts[grp1,42], col="red")
abline(a=0,b=theta_e@counts[grp2,7]/theta_e@counts[grp2,42], col="red")
plot(theta_e@counts[,42]/theta_e@counts[,7], col=ifelse(theta_e@group=="DMSO","red","blue"))

theta_f=setActive(pd, what = "theta_f")

### plot log-ratio abundance
parallel(theta_d, cutoff= .05, include="reference_Gene")
parallel(theta_e, cutoff= .05, include="reference_Gene")
###########################################################

