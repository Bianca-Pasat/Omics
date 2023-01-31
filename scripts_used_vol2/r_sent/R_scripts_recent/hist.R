boxplot(affyRaw,  target="core",las = 2, main = "Pre-Normalization")
affyRaw1 <- read.celfiles("/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_DN41.CEL","/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_DN42.CEL","/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_DN43.CEL")
affyRaw2 <- read.celfiles("/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_MN41.CEL","/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_MN42.CEL","/media/bianca/Elements/ASAG_Lab/home/bianca/CELL_FILES/MDA231_MN43.CEL" )
eset1<-rma(affyRaw1)
eset2<-rma(affyRaw2)
head(eset1)
eset1@assayData$exprs
exp<-cbind(eset1@assayData$exprs,eset2@assayData$exprs)
head(exp)
boxplot(Exp)
boxplot(exp)
write.table(exp,"/home/bianca/microarray_my_norm.txt", sep="\t", quote=F, row.names=rownames(exp))
source("~/.active-rstudio-document")
micbio<-read.table("/home/bianca/Downloads/microarray_UP_prioritized.txt", sep="\t")
rnabio<-read.table("/home/bianca/Downloads/rnaseq_UP_prioritized.txt", sep="\t")
ridd_all<-read.table("/home/bianca/Downloads/mrna-microarray-out_ids_all", sep="\t")
ridd_exact<-read.table("/home/bianca/Downloads/mrna-microarray-exact-seq", sep="\t")
# #annotation of genes
httr::set_config(httr::config(ssl_verifypeer = FALSE))
homo.anno = useEnsembl(biomart = "ensembl",
dataset = "hsapiens_gene_ensembl",
mirror = "useast")
library(biomaRt)
values=c(micbio$V2, rnabio$V2)
attributes <- c("ensembl_gene_id", "hgnc_symbol")
filters="ensembl_gene_id"
genes_ridd_prior <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
homo.anno = useEnsembl(biomart = "ensembl",
dataset = "hsapiens_gene_ensembl",
mirror = "useast")
genes_ridd_prior <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
genes_ridd_prior
filters
values
filters="hgnc_symbol"
genes_ridd_prior <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
genes_ridd_prior
head(ridd_all)
grp<-ridd_all[ridd_all$V1 %in% genes_ridd_prior$ensembl_gene_id,]
grp1<-merge(ridd_all, ridd_all, by=1)
dim(grp1)
dim(grp)
length(grp)
length(unique(grp))
length(values)
91/133
grp1<-merge(ridd_all, ridd_all, by.x=1, by.y=1)
head(grp1)
dim(grp1)
grp1<-merge(ridd_all, genes_ridd_prior, by.x=1, by.y=1)
dim(grp1)
dim(grp)
dim(unique(grp)0
dim(unique(grp))
length(unique(grp1))
head(grp``)
head(grp1)
dim(unique(grp1$V1))
length(unique(grp1$V1))
dim(unique(grp))
length(unique(grp))
grep<-ridd_exact[ridd_exact$V1 %in% genes_ridd_prior$ensembl_gene_id,]
grep1<-ridd_exact[ridd_exact$V1 %in% genes_ridd_prior$ensembl_gene_id,]
ugrp<-unique(grp)
grep1<-ridd_exact[ridd_exact$V1 %in% genes_ridd_prior$ensembl_gene_id,]
ugrep1<-unique(grep1)
dim(ugrep1)
length(ugrep1)
41/91
ugrep1
dim(grep1)
length(grep1)
grep1<-genes_ridd_prior[ridd_exact$V1 %in% genes_ridd_prior$ensembl_gene_id,]
length(grep1)
grep1
### bio info and exact and all ridd
grp<-merge(ridd_all, genes_ridd_prior, by.x=1, by.y=1)
head(grp)
grp<-merge(ridd_all, genes_ridd_prior, by.x=1, by.y=1)
ugrp<-unique(grp$hgnc_symbol)
grp1<-merge(ridd_exact, genes_ridd_prior, by.x=1, by.y=1)
ugrep1<-unique(grep1$ensembl_gene_id)
ugrp
BiocManager::install("ggVennDiagram")
library("ggVennDiagram")
a<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("up RNAseq","up Microarray")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
x1<-list(grp,ridd_all)
a<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("up RNAseq","up Microarray")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
x1<-list(as.character(grp$V1),as.character(ridd_all))
a<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("up RNAseq","up Microarray")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
a
head(ridd_a;;
head(ridd_all)
head(bio_common)
x1<-list(as.character(values),as.character(ridd_all$V1))
x2<-list(as.character(values),as.character(ridd_exact$V1))
a<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("Prioritized genes","Genes with CNGCNG")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
b<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("Prioritized genes","Genes with CUGCAG")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
a
values
x1<-list(as.character(genes_ridd_prior$ensembl_gene_id),as.character(ridd_all$V1))
x2<-list(as.character(genes_ridd_prior$ensembl_gene_id),as.character(ridd_exact$V1))
a<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("Prioritized genes","Genes with CNGCNG")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
a
b<-ggVennDiagram(
x1, label_alpha = 0,
category.names = c("Prioritized genes","Genes with CUGCAG")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
b
b<-ggVennDiagram(
x2, label_alpha = 0,
category.names = c("Prioritized genes","Genes with CUGCAG")
) +
ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")
b
rnaseq8<-read.table("/home/bianca/Downloads/rnaseq_rankprod_DEregulatedBio_8.txt", sep="\t")
rnaseq24<-read.table("/home/bianca/Downloads/rnaseq_rankprod_DEregulatedBio_24.txt", sep="\t")
library(biomaRt)
homo.anno = useEnsembl(biomart = "ensembl",
dataset = "hsapiens_gene_ensembl",
mirror = "useast")
valuesr8<-rownames(rnaseq8)
valuesr24<-rownames(rnaseq24)
attributes <- c("ensembl_gene_id", "hgnc_symbol")
filters="ensembl_gene_id"
genes_r8 <- getBM(attributes = attributes, filters=filters,values = valuesr8, mart = homo.anno)
genes_r24 <- getBM(attributes = attributes, filters=filters,values = valuesr24, mart = homo.anno)
genes_r24
head(genes_r8)
de=rnaseq8
de$ens<-rownames(de)
genesh<-genes_r8
de<-merge(geneseh, de, by.x=1, by.y="ens")
de<-merge(geneseh, de, by.x=1, by.y="ens")
geneseh<-genes_r8
de<-merge(geneseh, de, by.x=1, by.y="ens")
head(de)
library(ggrepel)
d0<-de
d0$diffexpressed <- "NO"
d0$diffexpressed[de$FC..class1.class2. > 0.5 & de$P.value < 0.05] <- "UP" ## BECAUSE RP DOES 0 VS 1
d0$diffexpressed[de$FC..class1.class2. < -0.5 & de$P.value < 0.05] <- "DOWN" ## BECAUSE RP DOES 0 VS 1
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed)) + geom_point() +
scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
#geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  + theme_minimal()
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed)) + geom_point() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4")) +
#geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  + theme_minimal()
ggplot(data=de, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed above 1")  +
theme_minimal()
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed above 1")  +
theme_minimal()
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed above 1")  +
theme_minimal()
# plot adding up all layers we have seen so far
ggplot(data=de, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=V2)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  +
theme_minimal()
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  +
theme_minimal()
genes_r24 <- getBM(attributes = attributes, filters=filters,values = valuesr24, mart = homo.anno)
head(genes_r24)
de=rnaseq24
de$ens<-rownames(de)
geneseh<-genes_r24
de<-merge(geneseh, de, by.x=1, by.y="ens")
d0<-de
d0$diffexpressed <- "NO"
d0$diffexpressed[de$FC..class1.class2. > 0.5 & de$P.value < 0.05] <- "UP" ## BECAUSE RP DOES 0 VS 1
d0$diffexpressed[de$FC..class1.class2. < -0.5 & de$P.value < 0.05] <- "DOWN" ## BECAUSE RP DOES 0 VS 1
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  +
theme_minimal()
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
geom_point() +
theme_minimal() +
geom_text_repel() +
scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  +
geom_vline(xintercept=c(-0.5, 0.5), col="red") +
geom_hline(yintercept=-log10(0.05), col="red") +
labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
ggtitle("Significantly differentially expressed with FC above 0.5")  +
theme_minimal()
load("/home/bianca/Downloads/exons.RData")
library(sva)
library(tidyverse)
library("DEXSeq")
library("reshape2")
library("ggplot2")
library("stringr")
head(counts)
head(a)
co<-dxd1@modelFrameBM$count
head(co)
co<-dxd1@assays@data
head(co)
co<-dxd1@assays@data$counts
head(co)
head(counts)
head(design\)
design
sampleTable
head(a)
co<-counts
counts<-a
counts<-a[-c(1:5),]
#colnames(counts)<-counts[length(counts[[1]]),]
#counts<-counts[-length(counts[[1]]),]
rownames(counts)<-counts$V1
counts=counts[,-1]
colnames(counts)<-sampleTable$sample
counts[,1:6]<-lapply(counts[1:6],as.numeric)
# function.norm <- function(counts){
#   norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE)
#   rownames(norm)=rownames(counts)
#   colnames(norm)=colnames(counts)
#
#   return(norm)
# }
# norm_samples = function.norm(counts)
counts1<-counts
keep.exprs <- aveLogCPM(counts) > 0
counts <- counts[keep.exprs,]
trans_cts<-as.data.frame(counts1)
# load packages
library(tidyverse)
sampleTable
sample_info<-sampleTable
pca_matrix <- trans_cts %>%
as.matrix() %>%
t()
# Perform the PCA
sample_pca <- prcomp(pca_matrix)
# Convert matrix to tibble
tibpca<-as_tibble(pca_matrix)
# Convert matrix to tibble - add colnames to a new column called "gene"
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
# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x
pc_scores <- pc_scores %>%
as_tibble(rownames = "sample")
# pc_scores %>%
#   # create the plot
#   ggplot(aes(PC1, PC2, col=sample_info$condition)) +
#   geom_point(size=3) +
#   geom_text(aes(label = sample), vjust = -1, nudge_y = 1) +
#   labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) + # substitute with % variances
#   ggtitle("PCA. by replicate") +
#   theme_bw()
# plot together
a<-sample_pca$x %>%
# convert it to a tibble
as_tibble(rownames = "sample") %>%
# join with "sample_info" table
full_join(sample_info, by = "sample") %>%
# make the plot
ggplot(aes(x =PC1, y = PC2, col = treatment, shape = replicate)) +
scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
geom_point(size=5) +
geom_text(aes(label = sample), vjust = -1, nudge_y = 1) +
labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
ggtitle("PCA of treatment and replicate") +
theme_bw()
a
sampleTable
a<-sample_pca$x %>%
# convert it to a tibble
as_tibble(rownames = "sample") %>%
# join with "sample_info" table
full_join(sample_info, by = "sample") %>%
# make the plot
ggplot(aes(x =PC1, y = PC2, col = condition, shape = replicate)) +
scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
geom_point(size=5) +
geom_text(aes(label = sample), vjust = -1, nudge_y = 1) +
labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
ggtitle("PCA of treatment and replicate") +
theme_bw()
a
trans_cts<- as.data.frame(counts1)
sample_info<-sampleTable
sample_info
raw_cts_long <- trans_cts %>%
pivot_longer(1:6, names_to = "sample", values_to = "expression")
# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
raw_cts_long <- trans_cts %>%
pivot_longer(DMSO_8_1:MKC8866_24_1, names_to = "sample", values_to = "expression")
# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
sample_info
raw_cts_long %>%
ggplot(aes(as.character(condition), log2(expression))) +
geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) +
facet_grid(cols = vars(replicate)) +
labs (y = "log2 expression", x = "time points") +
ggtitle("Boxplot of raw counts before normalization") +
theme_bw()
sample_info
raw_cts_long <- trans_cts %>%
pivot_longer(1:6, names_to = "sample", values_to = "expression")
# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
raw_cts_long
raw_cts_long %>%
ggplot(aes(as.character(condition), log2(expression))) +
geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) +
facet_grid(cols = vars(replicate)) +
labs (y = "log2 expression", x = "time points") +
ggtitle("Boxplot of raw counts before normalization") +
theme_bw()
raw_cts_long %>%
ggplot(aes(as.character(condition), log2(expression),fill = treatment)) +
geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) +
facet_grid(cols = vars(replicate)) +
labs (y = "log2 expression", x = "time points") +
ggtitle("Boxplot of raw counts before normalization") +
theme_bw()
raw_cts_long %>%
ggplot(aes(as.character(condition), log2(expression),fill = condition)) +
geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) +
facet_grid(cols = vars(replicate)) +
labs (y = "log2 expression", x = "time points") +
ggtitle("Boxplot of raw counts before normalization") +
theme_bw()
keep.exprs <- aveLogCPM(counts) > 0
counts <- counts[keep.exprs,]
library(edgeR)
keep.exprs <- aveLogCPM(counts) > 0
counts <- counts[keep.exprs,]
trans_cts<- as.data.frame(counts)
raw_cts_long <- trans_cts %>%
pivot_longer(1:6, names_to = "sample", values_to = "expression")
# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
raw_cts_long %>%
ggplot(aes(as.character(condition), log2(expression),fill = condition)) +
geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) +
facet_grid(cols = vars(replicate)) +
labs (y = "log2 expression", x = "time points") +
ggtitle("Boxplot of raw counts after normalization") +
theme_bw()
### VISUALIZATION
# MA plot
plotMA( dxr1, cex=0.8 )
dxr1
length(which(dxr1$padj < 0.01))
length(which(dxr1$padj < 0.01 && dxr1$log2fold_MKC_DMSO > 1))
length(which(dxr1$padj < 0.01 && dxr1$log2fold_MKC_DMSO > 0.5))
length(which(dxr1$padj < 0.01 && dxr1$log2fold_MKC_DMSO > 0))
length(which(dxr1$padj < 0.01 && which(dxr1$log2fold_MKC_DMSO > 0)))
now<-dxr1[dxr1$padj < 0.01,]
head(hitsire1WOW)
View(ire1_exon_seq)
### VISUALIZATION
# MA plot
plotMA( dxr1, cex=0.8 )
mongo<-read.table("/home/bianca/Downloads/interactions.tsv",sep="\t",quote=F)
mongo<-read.table("/home/bianca/Downloads/interactions.tsv",sep="\t")
mongo<-read.table("/home/bianca/Downloads/interactions.tsv",sep=" ")
mongo<-read.table("/home/bianca/Downloads/interactions.tsv",sep=",")
mongo<-read.table("/home/bianca/Downloads/interactions.tsv")
library(data.table)
mongo<-fread("/home/bianca/Downloads/interactions.tsv",data.table=F
)
head(mongo)
library(tidyverse)
mongonow<-mongo %>% select(c(1,8)) %>% group_by(2)
head(mongonow)
mongonow<-mongo %>% select(c(1,8)) %>% group_by(mongo$gene_name)
head(mongonow)
dim(mongonow)
dim(mongo)
write.table(mongo[,c(1,9),"/home/bianca/Downloads/interactions_small.csv",sep=",",quote=F,row.names=F]
)
write.table(mongo[,c(1,9)],"/home/bianca/Downloads/interactions_small.csv",sep=",",quote=F,row.names=F)
colnames(mongo)
write.table(mongo[,c(1,8)],"/home/bianca/Downloads/interactions_small.csv",sep=",",quote=F,row.names=F)
nala<-as.data.frame(mongo[,c(1,8)])
dimn(nala)
dim(nala)
head(nala)
write.table(nala,"/home/bianca/Downloads/interactions_small.csv",sep=",",quote=F,row.names=F)
?write.table
nala2<-nala %>% pivot_wider(,names_from = gene_name)
nala2<-nala %>% pivot_wider(,names_from = "gene_name")
nala2<-nala %>% pivot_wider(names_from =gene_name)
nala2<-nala %>% pivot_wider(names_from = nala$gene_name)
nala2<-nala %>% pivot_wider(names_to=drug, values_from = drug_name )
nala2<-nala %>% pivot_wider(names_from = gene_name, values_from = drug_name )
?pivot_longer
mongo<-fread("/home/bianca/Downloads/interactions.tsv",data.table=F,na.strings = " ")
which(is.na(mongo$drug_name))
mongo<-fread("/home/bianca/Downloads/interactions.tsv",data.table=F,na.strings = "")
which(is.na(mongo$drug_name))
nala<-as.data.frame(mongo[,c(1,8)])
nala2<-nala %>% pivot_wider(names_from = gene_name, values_from = drug_name )
head(nala2)
nala2
dim(nala2)
pivot_wider()
?pivot_wider()
length(unique(mongo$gene_name))
nala2<-nala %>% pivot_wider(names_from = drug_name, values_from = gene_name )
dim(nala2)
nala2<-nala %>% group_by_at(1))
nala2<-nala %>% group_by_at(1)
dim(nala2)
head(nala)
dim(nala)==dim(nala2)
nala2<-nala %>% pivot_wider(names_from = gene_name, values_from = drug_name )
nala2[8,]
write.table(nala2,"/home/bianca/Downloads/tible_interactions.csv",sep=",",quote=F, row.names = F)
df<-apply(nala2,2,as.character)
hewad(df)
head(df)
write.csv(df,"/home/bianca/Downloads/tible_interactions.csv",sep=",",quote=F, row.names = F)
write.csv(df,"/home/bianca/Downloads/tible_interactions.csv",quote=F, row.names = F)
df<-apply(nala2,c(1,2),as.character)
head(df)
write.csv(df,"/home/bianca/Downloads/tible_interactions.csv",quote=F, row.names = F)
head(nala2)
nala2[[1]]
savehistory("~/hist.R")
