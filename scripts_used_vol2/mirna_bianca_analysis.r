library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library(stringr)

# 
# raw_counts <- read.delim("/home/bianca/Downloads/reads_one_per_mirna.csv", sep=",", quote = "") 
# raw_counts <- read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/miRNA-results/reads_one_per_mirna.csv", sep=",", quote = "") 
# precursor<-raw_counts[,4]
# 
# rownames(raw_counts)<-raw_counts[,3]
# exp <- raw_counts[,c(5:10,17:22)]
# colnames(exp)<-str_replace(colnames(exp), "[X][.]", "")
# colnames(exp)<-str_replace(colnames(exp), "[.]", "")
# rownames(exp)<-str_replace(rownames(exp), "[\"]", "")
# rownames(exp)<-str_replace(rownames(exp), "[\"]", "")

## START HERE
exp <- read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/miRNA-results/mirnas.csv", sep=",", quote = "", row.names = 1, header=TRUE) 
rox<-rownames(exp)
X=as.data.frame(lapply(exp, function(x) as.numeric(x)))
rownames(X)=rownames(exp)
exp<-X

sample_info<-data.frame(colnames(exp))
sample_info$condition<-c(rep(c(rep(0,3), rep(1,3)),2))
sample_info$replicate<-as.factor(c(rep(c("1","2","3"), 4)))
sample_info$times<-c(rep(8,6), rep(24,6))
sample_info$treatmentL<-c(rep(paste("DMSO","8",sep="_"),3),rep(paste("MKC","8",sep="_"),3),rep(paste("DMSO","24",sep="_"),3),rep(paste("MKC","24",sep="_"),3))
sample_info$treatment<-as.factor(c(rep(c(rep("DMSO",3), rep("MKC",3)),2)))
colnames(sample_info)<-c("sample", "condition", "replicate", "times","treatment","drug")

y <- DGEList(counts=exp, samples=sample_info , genes=rownames(exp))

## com and log2 
y$counts <- cpm(y, log=TRUE)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- rowSums(y$counts==0) < 7
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

### or

# FILTRATION OF LOW COUNTS
# row has at least two columns with a count of over 2
# keep=apply(y$counts, 1, function(x) { length(x[x>2])>=7 } )
# Y= y$counts[keep,]

### QUANTILE NORMALIZATION (or whichever)
#normalization


function.norm <- function(counts){
  norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
  rownames(norm)=rownames(counts)
  colnames(norm)=colnames(counts)
  
  return(norm)
}

norm_samplesm = function.norm(y$counts) # this one

norm_samplesm = function.norm(log2(y$counts))
norm_samples2 = function.norm(log2(Y))
norm_samples = function.norm(as.matrix(y$counts)) # no log)

################################################
# batch exploration
## on all
norm_samplesm=norm_samples2
PreparePlotpca(norm_samplesm, sample_info,sample_info$treatment,as.factor(sample_info$times),values=cbp_12)
PreparePlotpca(norm_samplesm, sample_info,as.factor(sample_info$treatment),as.factor(sample_info$replicate),values=cbp_12)

comsvaRe<-comsva(norm_samplesm, sample_info, treatment, replicate)
PreparePlotpca(comsvaRe, sample_info,sample_info$treatment,as.factor(sample_info$times),values=cbp_12)
PreparePlotpca(comsvaRe, sample_info,sample_info$condition,as.factor(sample_info$times),values=cbp_12)

svacomre<-svacom(as.matrix(norm_samplesm),sample_info,treatment,replicate)
PreparePlotpca(svacomre, sample_info,sample_info$treatment,as.factor(sample_info$times),values=cbp_12)
PreparePlotpca(svacomre, sample_info,sample_info$condition,as.factor(sample_info$times),values=cbp_12)

comre<-com(norm_samplesm,sample_info,replicate)
PreparePlotpca(comre, sample_info,sample_info$treatment,as.factor(sample_info$times),values=cbp_12)
PreparePlotpca(comre, sample_info,sample_info$condition,as.factor(sample_info$times),values=cbp_12)

savre<-svAll(as.matrix(norm_samplesm),sample_info,treatment,replicate)
PreparePlotpca(savre, sample_info,sample_info$treatment,as.factor(sample_info$times),values=cbp_12)
PreparePlotpca(savre, sample_info,sample_info$condition,as.factor(sample_info$times),values=cbp_12)


clean_mirna_sva<-savre
clean_mirna_comsva<-comsvaRe

saveRDS(clean_mirna_sva,"/home/ben/Desktop/mirna_bianca/clean_mirna_sva")
saveRDS(clean_mirna_comsva,"/home/ben/Desktop/mirna_bianca/clean_mirna_comsva")


### DIFFERENTIAL EXPRESSION 
### # RP parameters:

is_logged <- as.logical("TRUE")
logbase <- 2
is_normalized <- as.logical("TRUE")
data_size=as.logical("FALSE")


#8 HOURS
samples<-clean_mirna_comsva[,c(1:6)]
cl=c(rep(1,3),rep(0,3))
out_dir="/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8.txt"
upsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_upregulatedBio.txt"
downsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_downregulatedBio.txt"
allsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt"

#24 HOURS
samples<-clean_mirna_comsva[,c(7:12)]
cl=c(rep(1,3),rep(0,3))
out_dir="/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_24.txt"
upsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_24_upregulatedBio.txt"
downsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_24_downregulatedBio.txt"
allsave= "/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_24_DEregulatedBio.txt"


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


#### LOAD DE AND SUBSTRACT FROM EACH TIME POINT!
de8<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8.txt",sep="\t")
de8<-rownames(de8)
de24<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_24.txt",sep="\t")
de24<-rownames(de24)

deAll<-c(de8,de24)

cleand_de_mirna<-clean_mirna_comsva[rownames(clean_mirna_comsva) %in% deAll,]
cleand_de_mirna8<-cleand_de_mirna[,c(1:6)]
cleand_de_mirna24<-cleand_de_mirna[,c(7:12)]

saveRDS(cleand_de_mirna,"/home/ben/Desktop/mirna_bianca/cleaned_de_mirna")
saveRDS(cleand_de_mirna8,"/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_8")
saveRDS(cleand_de_mirna24,"/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_24")

filnorm_mirna8<-norm_samplesm[rownames(norm_samplesm) %in% de8,]
filnorm_mirna8<-filnorm_mirna8[,c(1:6)]
filnorm_mirna24<-norm_samplesm[rownames(norm_samplesm) %in% de24,]
filnorm_mirna24<-filnorm_mirna24[,c(7:12)]
saveRDS(filnorm_mirna8,"/home/ben/Desktop/mirna_bianca/filnorm_mirna_8")
saveRDS(filnorm_mirna24,"/home/ben/Desktop/mirna_bianca/filnorm_mirna_24")


# #Load the mapping info
# load('hsa-vtm-gene.Rdata')
# #the mapping can be obtained using SpidermiR package
# 
# #Identify the miRBase version of mapping file and data
# library(miRNAmeConverter)
# library(SpidermiR)
# 
# #miRNAs_mapping<-SpidermiRdownload_miRNAprediction(mirna_list=deAll)
# saveRDS(miRNAs_mapping,"/home/ben/Desktop/mirna_bianca/mapping")
# miRNAs_mapping=readRDS("/home/ben/Desktop/mirna_bianca/mapping")
# miRNAs <- deAll
# 
# nc <- MiRNANameConverter()
# 
# # assessVersion(nc, miRNAs, verbose = FALSE) #version 9.2
# # assessVersion(nc, miRNAs_mapping, verbose = FALSE) #version 20
# 
# #convert miRNA annotation from our miRNA matrix from v.9.2 to v.20
# # library(anamiR)
# # new_mirna_matrix <- miR_converter(mirna_gbm_corrected, remove_old = TRUE, original_version=9.2, latest_version = 20)
# 
# #restricting mapping info to measured miRNA
# idx<-intersect(miRNAs_mapping$V1, miRNAs)
# 
# id_mapping <- miRNAs_mapping[miRNAs_mapping$V1 %in% miRNAs,] 
# #creating the miRNA association tables
# mirna2GeneMap <- stack(deAll) #convert in a data.frame. Each line has a gene and one mirna. 
# 
# mirna2GeneMap=id_mapping 
# names(mirna2GeneMap) <- c('miRNA','Gene')
# head(mirna2GeneMap)
# #creating the mRNA association tables
# 
# de_rna_8<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")
# de_rna_24<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_24")
# 
# de_rna_all<-c(rownames(de_rna_8), rownames(de_rna_24))
# 
# expr2GeneMap <- data.frame(measurement = de_rna_all, 
#                            Gene = de_rna_all)
# 
# expr2GeneMap <- data.frame(measurement = rownames(rna_gbm_corrected), 
#                            Gene = rownames(rna_gbm_corrected))
# 
# head(expr2GeneMap)
# #creating the data mapping
# dataMappingExprMirna <- combiningMappings(mappings = list(expr = expr2GeneMap, 
#                                                           mirna = mirna2GeneMap),
#                                           retainAll = TRUE, reference = 'Gene'); 
# 
# dataMappingExprMirna2 <- combiningMappings(mappings = list(expr = expr2GeneMap, 
#                                                           mirna = mirna2GeneMap),
#                                           retainAll = FALSE, reference = 'Gene'); 
# head(dataMappingExprMirna)
# 
# #Overlaping samples
# #restrict the datasets to the elements of the data mapping and those samples that are common between omics
# 
# # 8 hours clean matrices
# sample8<-sample_info[1:6,]
# rownames(sample8)=sample8$sample
# 
# exprTMP=readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")
# mirnaTMP=readRDS("/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_8")
# 
# #specifying the data types.
# 
# dataTypesExprMirna <- c("count", "count")
# 
# 
# #preparing the datasets
# exprTMP <- createOmicsExpressionSet(Data = exprTMP)
# 
# 
# mirnaTMP <- createOmicsExpressionSet(Data = as.matrix(mirnaTMP),sample8)
# 
# 
# dataInputExprMirna <- list(expr = exprTMP, mirna = mirnaTMP)
# #Parametric Combination
# set.seed(12345)
# omicsPCRes_cPC <- omicsNPC(dataInput = dataInputExprMirna,
#                           dataMapping = dataMappingExprMirna,
#                           dataTypes = dataTypesExprMirna,
#                           verbose = TRUE)

