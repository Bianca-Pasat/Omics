library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library(stringr)


raw_counts <- read.delim("/home/bianca/Downloads/reads_one_per_mirna.csv", sep=",", quote = "") 
raw_counts <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/miRNA-results/reads_one_per_mirna.csv", sep=",", quote = "") 
precursor<-raw_counts[,4]
rownames(raw_counts)<-raw_counts[,3]
exp <- raw_counts[,c(5:10,17:22)]
colnames(exp)<-str_replace(colnames(exp), "[X][.]", "")
colnames(exp)<-str_replace(colnames(exp), "[.]", "")
rownames(exp)<-str_replace(rownames(exp), "[\"]", "")
rownames(exp)<-str_replace(rownames(exp), "[\"]", "")
sample_info<-data.frame(colnames(exp))
sample_info$condition<-c(rep(c(rep(0,3), rep(1,3)),2))
sample_info$replicate<-as.factor(c(rep(c("1","2","3"), 4)))
sample_info$times<-c(rep(8,6), rep(24,6))
sample_info$treatmentL<-c(rep(paste("DMSO","8",sep="_"),3),rep(paste("MKC","8",sep="_"),3),rep(paste("DMSO","24",sep="_"),3),rep(paste("MKC","24",sep="_"),3))
sample_info$treatment<-as.factor(c(rep(c(rep("DMSO",3), rep("MKC",3)),2)))
colnames(sample_info)<-c("sample", "condition", "replicate", "times","treatment","drug")
samples_info_mirna<-sample_info


#mirna<-read.table("/home/bianca/mirna_aitor.txt", sep="\t")
#exp<-exp1

y <- DGEList(counts=exp, samples=sample_info , genes=rownames(exp))


# REMOVE LOW EXPRESSED GENES
keep.exprs <- rowSums(y$counts==0) < 7
y <- y[keep.exprs,, keep.lib.sizes=FALSE]


#normalization

function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  
  return(norm)
}


norm_samplesm = function.norm(log2(y$counts))
norm_samples = function.norm(y$counts)

################################################
# batch exploration
## on all
PreparePlotpca(mstus_lipids, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(mstus_lipids, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

comsvaRe<-comsva(mstus_lipids, sampleInfo, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(comsvaRe, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

svacomre<-svacom(as.matrix(mstus_lipids),sampleInfo,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(svacomre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

comre<-com(mstus_lipids,sampleInfo,replicate)
PreparePlotpca(comre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(comre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

savre<-svAll(as.matrix(mstus_lipids),sampleInfo,treatment,replicate)
PreparePlotpca(savre, sampleInfo,sampleInfo$treatment,as.factor(sampleInfo$Time_point),values=cbp_12)
PreparePlotpca(savre, sampleInfo,sampleInfo$condition,as.factor(sampleInfo$Time_point),values=cbp_12)

## time points

# 6 hours
sampleInfo6<-sampleInfo[which(sampleInfo$Time_point =="6 hour"),]
samples6<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo6$sample]

PreparePlotpca(samples6, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

comsvaRe<-comsva(samples6, sampleInfo6, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo6,sampleInfo6$condition,as.factor(sampleInfo6$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples6),sampleInfo6,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

comre<-com(samples6,sampleInfo6,replicate)
PreparePlotpca(comre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples6),sampleInfo6,treatment,replicate)
PreparePlotpca(savre, sampleInfo6,sampleInfo6$treatment,as.factor(sampleInfo6$replicate),values=cbp_12)

clean6<-comsvaRe
saveRDS(clean6, "/home/ben/Desktop/lipidsb/clean6")

## 12
sampleInfo12<-sampleInfo[which(sampleInfo$Time_point =="12 hour"),]
samples12<-mstus_lipids[,colnames(mstus_lipids) %in% sampleInfo12$sample]

PreparePlotpca(samples12, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

comsvaRe<-comsva(samples12, sampleInfo12, treatment, replicate)
PreparePlotpca(comsvaRe, sampleInfo12,sampleInfo12$condition,as.factor(sampleInfo12$replicate),values=cbp_12)

svacomre<-svacom(as.matrix(samples12),sampleInfo12,treatment,replicate)
PreparePlotpca(svacomre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

comre<-com(samples12,sampleInfo12,replicate)
PreparePlotpca(comre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)

savre<-svAll(as.matrix(samples12),sampleInfo12,treatment,replicate)
PreparePlotpca(savre, sampleInfo12,sampleInfo12$treatment,as.factor(sampleInfo12$replicate),values=cbp_12)


# 8 hours
exp1<-y$counts[,c(1:6)]
samples_info<-sample_info[c(1:6),]
y <- DGEList(counts=exp1, samples=samples_info , genes=rownames(exp1))
sample_info<-samples_info


# 24 hours
exp2<-y$counts[,c(6:12)]
samples_info<-sample_info[c(6:12),]
y <- DGEList(counts=exp2, samples=samples_info , genes=rownames(exp2))
sample_info<-samples_info



### OTHER NORMALIZATION METHODS

# # COUNT PER MILLION, LOG2 COUNTS. for plotting 
# lcpm_pre <- cpm(y, log=TRUE) # keep this to compare the normalization method afterwards
# 
# y$counts <- cpm(y, log=F) 
# #Normalize 
# y <- calcNormFactors(y)




