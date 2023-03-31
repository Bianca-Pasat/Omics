

### DEXSeq Differential Exon Usage Analysis
# initial steps of analysis are done using two Python scripts provided with DEXSeq
#
# pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
# list.files(pythonScriptsDir)

# system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

###  Preparing the annotation
# ubuntu command line shell
# Make sure that your current working directory contains the GTF file 
# python path_to_script GTF_file_name output_file.gff
# ex. python /mnt/c/Users/alman/Documents/R/win-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py Drosophila_melanogaster.BDGP5.72.gtf Dmel.BDGP5.25.62.DEXSeq.chr.gff

###  Counting reads per exon
# for each sam/bam file
# python /path/to/dexseq_count.py output_file.gff untreated1.sam untreated1fb.txt 
# -f bam
# -s reverse

### Reading the data in to R 
library("DEXSeq")
library("reshape2")
library("ggplot2")
library("stringr")

setwd("/home/ben/Desktop/EXONS_24_AITOR/")
countFiles = list.files( pattern=".txt$", full.names=TRUE)
countFiles <- countFiles[c(1:3,7:9)]
flattenedFile = list.files( pattern="gtf", full.names=TRUE)

sampleTable = data.frame(
  row.names = c( "D24R1", "D24R2", "D24R3", 
                 "M24R1", "M24R2", "M24R3" ),
  condition = c("DMSO", "DMSO", "DMSO",  
                "MKC", "MKC", "MKC" ),
  libType = c( "single-end", "single-end", "single-end", 
               "single-end", "single-end", "single-end" ))

#construct an DEXSeqDataSet object from this data
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon)

# first 6 columns = number of reads mapping to exonic regions, last 6 = sum of the counts for rest of the exons from the same gene on each sample.

### Normalization (same method as in DEseq2)
dxd = estimateSizeFactors( dxd )

# Dispersion estimation (estimate strength of the noise)
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )

### Testing for differential exon usage
# fits a generalized linear model with the formula ~sample + exon + condition:exon and compare it to the smaller model (the null model) ~ sample + exon.

dxd = testForDEU( dxd )

# estimate relative exon usage fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

### Results summary
dxr1 = DEXSeqResults( dxd )
mcols(dxr1)$description
resultados <- as.data.frame(dxr1)

DEXSeqHTML(dxr1)


##### LOAD MATRICES THAT DON'T CONTAIN 0: REMOVED LOWLY EXPRESSED VALUES
library(omicade4)
bibi<-function(x_matrix, pars){
  x2<-data.frame(row.names = rownames(x_matrix))
  for (i in 1:length(pars$hours)){
    col1=which(x_info$hours==pars$hours[i])
    col2=which(x_info$treatment=="DMSO")
    col3=which(x_info$treatment!="DMSO")
    x1<-x_matrix %>% mutate(meanBase=rowMeans(x_matrix[,col2]))
    x1<-sapply(x1[,col3], function(x) {x - x1$meanBase})
    x2<-cbind(x2,x1)
    return(x2)
  }}

# mrna
mrna8<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rna_filt_norm_de_8"))
mrna24<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rna_filt_norm_de_24"))
demrna<-c(rownames(mrna8), rownames(mrna24))
fnmrna<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mrna/rna_normalized_final.txt")
x_matrix<-fnmrna[demrna,]
write.table(x_matrix,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesfilterednorm_mrna.txt", quote=F, row.names=T, sep = "\t")

x_info<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mrna/all_sample_info.txt", sep="\t")
colnames(x_info)<-x_info[1,]
x_info<-x_info[-1,]
x_matrix<-x_matrix[,x_info$sample]
pars<-list(hours=c(8,24),
           treatment="DMSO")
mrna<-bibi(x_matrix, pars)
write.table(mrna,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesMKC-DMSO_mrna.txt", quote=F, row.names=T, sep = "\t")

## proteins
deprot<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_all_de_prot.txt")
fnprot<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/proteins/removed_low_counts_proteins_QUANTILE_NORM_log2"))
x_matrix<-fnprot[deprot$V1,]
write.table(x_matrix,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesfilterednorm_prot.txt", quote=F, row.names=T, sep = "\t")

x_info<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/proteins/proteins_sample_info"))
x_info$condition<-x_info$treatment
x_info$treatment<-gsub("_.*","",x_info$treatment)
colnames(x_info)[3]<-"hours"
x_matrix<-x_matrix[,x_info$sample]

pars<-list(hours=c(8,24,48,72),
           treatment="DMSO")
prot<-bibi(x_matrix, pars)
write.table(prot,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesMKC-DMSO_prot.txt", quote=F, row.names=T, sep = "\t")

# mirna 
demirna<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/mirna/cleaned_de_mirna"))
fnmirna<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/filtered_norm_all.txt")
x_matrix<-fnmirna[rownames(demirna),]
write.table(x_matrix,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesfilterednorm_mirna.txt", quote=F, row.names=T, sep = "\t")

x_info<-readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/mirna/sample_info")
x_matrix<-x_matrix[,x_info$sample]

pars<-list(hours=c(8,24),
           treatment="DMSO")
mirna<-bibi(x_matrix, pars)23+
write.table(mirna,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesMKC-DMSO_mirna.txt", quote=F, row.names=T, sep = "\t")

## lipids
delips<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod_all_de_lipids.csv")
fnlipid<-as.data.frame(readRDS("/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/lipids/mstus_lipids_all"))
x_matrix<-fnlipid[delips$V1,]
write.table(fnlipid,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesfilterednorm_lipid.txt", quote=F, row.names=T, sep = "\t")

x_info<-as.data.frame(readRDS("//home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/lipids/lipids_sample_info"))
# colnames(x_info)<-c("sample","hours","treatment","replicate")
# x_info$hours<-gsub(" .*","",x_info$hours)
# saveRDS(x_info,"/home/ben/Desktop/OMICS_JANUARY_2023/inputs/filtered_norm/lipids/lipids_sample_info")
x_matrix<-x_matrix[,x_info$sample]

pars<-list(hours=c(6,12,24,48,72),
           treatment="DMSO")
lipids<-bibi(x_matrix, pars)
write.table(lipids,"/home/ben/Desktop/OMICS_JANUARY_2023/rankprodDE_valuesMKC-DMSO_prot.txt", quote=F, row.names=T, sep = "\t")
# LIPIDS NEED MORE WORK TO FIGURE WHICH SAMPLES AND HOW

#####

prot<-readRDS("/home/bianca/Desktop/prot.RDS")
x_b=rownames(prot)
prot1=as.data.frame(sapply(prot, function(x) as.numeric(x)))
rownames(prot)=x_b

#rna<-readRDS("/home/bianca/Desktop/rna.RDS")
rna<-readRDS("/home/bianca/Desktop/rna_filt_norm.RDS")
#rna=rna[,c(1,5,4,6,3,2,9,8,10,12,11,7)]
#rna<-read.table("/home/bianca/Downloads/rna_normalized_final.txt", sep="\t")
#rna=rna[,c(1,5,4,6,3,2,9,8,10,12,11,7)]
#saveRDS(rna,"/home/bianca/Desktop/rna_filt_norm.RDS")
#write.table(sample_info_RNAseq,"/home/bianca/Downloads/rna_normalized_final_sample_info.txt", sep="\t",quote=F,  row.names=T)

#### FIX MICROARRAY TO HAVE HGNC SYMBOL!!!!!!!!
micro<-readRDS("/home/bianca/Desktop/micro.RDS")
#colnames(micro)=rep("MKC_24",3)


pr1=prot[,c(1,2,3,7,8,9)]
r1=rna[,c(1,2,3,7,8,9)]
mda<-list(as.matrix(pr1), as.matrix(r1))
mda<-list(as.matrix(prot), as.matrix(rna))

####
mda<-list(as.matrix(prot), as.matrix(mrna))

#### all 8 de
demirna8<-read.table("/home/ben/Desktop/mirna_bianca/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
colnames(demirna8)<-c("log-ratio","Pvalue")
demirna8$Gene=rownames(demirna8)
demirna8=demirna8[,c(3,1,2)]

mirna8<-readRDS("/home/ben/Desktop/mirna_bianca/cleaned_de_mirna_8")
mrna8<-readRDS("/home/ben/Desktop/rna_bianca/clean_DeRna_8")

demi8<-mirna8[rownames(mirna8) %in% demirna8$Gene, ]
dedemi8<-demirna8[demirna8$Gene %in% rownames(mirna8) , ]

mrna8<-mrna8[,colnames(mirna8)]

prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_de_prot8")
lipids6<-readRDS("/home/ben/Desktop/lipidsb/lipids6")
lipred<-lipids6[,c(1,3,4,5,6,8)]
colnames(lipred)=colnames(prot8)


### start with clean datasets but not de
mrna8=readRDS("/home/ben/Desktop/rna_bianca/clean_AllRna_8")
mirna8=readRDS("/home/ben/Desktop/mirna_bianca/clean_mirna_comsva")
mirna8<-mirna8[,c(1:6)]
mrna8<-mrna8[,colnames(mirna8)]
lipids6<-readRDS("/home/ben/Desktop/lipidsb/clean6")
lipred<-lipids6[,c(1,3,4,5,6,8)]
prot8<-readRDS("/home/ben/Desktop/proteomics_cleaned_bianca/clean_prot8")
colnames(lipred)=colnames(prot8)


mda<-list(as.matrix(mrna8), as.matrix(demi8))
mda<-list(mRNA=as.matrix(mrna8), miRNA=as.matrix(demi8), lipids=as.matrix(lipred), proteins=as.matrix(prot8))

mda<-list(mRNA=as.matrix(as.data.frame(mrna)),  proteins=as.matrix(as.data.frame(prot)), miRNA=as.matrix(as.data.frame(mirna)))



all(apply((x <- sapply(mda, colnames)), 2, function(y)
  + identical(y, x[,1])))


layout(matrix(1:4, 2, 2))
par(mar=c(2, 1, 0.1, 6))
for (df in mda) {
  d <- dist(t(df))
  hcl <- hclust(d)
  dend <- as.dendrogram(hcl)
  plot(dend, horiz=TRUE)
}


mcoin <- mcia(mda)
cond=c(rep("DMSO",3),rep("MKC",3))
### plot to see which genes differentiate your conditions!
plot(mcoin, axes=1:2, phenovec=cond, sample.lab=FALSE, df.color=1:2)
plot(mcoin, axes=1:2, phenovec=cond, sample.lab=FALSE, df.color=1:4)

## choose genes~proteins according to your pca
goi <- selectVar(mcoin, a1.lim=c(-Inf,-2), a2.lim=c(-Inf, Inf))
goi <- selectVar(mcoin, a1.lim=c(-Inf,), a2.lim=c(-Inf, Inf))

goi2 <- selectVar(mcoin, a1.lim=c(2,Inf), a2.lim=c(-Inf, Inf))
goi05 <- selectVar(mcoin, a1.lim=c(0.5,Inf), a2.lim=c(-Inf, Inf))
goi1 <- selectVar(mcoin, a1.lim=c(1,Inf), a2.lim=c(-Inf, Inf))
goi0 <- selectVar(mcoin, a1.lim=c(0,Inf), a2.lim=c(-Inf, Inf))

write.table(goi2$var,"/home/ben/Desktop/integrated_mcia_2_features", sep="\t", quote=F, row.names = F)
write.table(goi05$var,"/home/ben/Desktop/integrated_mcia_05_features", sep="\t", quote=F, row.names = F)
write.table(goi1$var,"/home/ben/Desktop/integrated_mcia_1_features_with_all_features", sep="\t", quote=F, row.names = F)
write.table(goi0$var,"/home/ben/Desktop/integrated_mcia_0_features", sep="\t", quote=F, row.names = F)


geneStat <- plotVar(mcoin, var=c("S100B", "S100A1"), var.lab=TRUE)
