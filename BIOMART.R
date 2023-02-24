if (!require("BiocManager", quietly = TRUE))
  install.packages("tidyverse")
BiocManager::install(version = "3.14")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("")

library(biomaRt)

frp4<-read.table("/home/bianca/Downloads/rankprod_results_4_for.txt", sep="\t")


rnaseq8<-read.table("/home/bianca/Downloads/rnaseq_rankprod_DEregulatedBio_8.txt", sep="\t")
rnaseq24<-read.table("/home/bianca/Downloads/rnaseq_rankprod_DEregulatedBio_24.txt", sep="\t")

mart = useDataset("hsapiens_gene_ensembl", mart=useMart("ensembl"))

listAttributes(mart)

homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")

ncbi_to_gene<-getBM(attributes = c("affy_huex_1_0_st_v2","hgnc_symbol"), 
                    filters="affy_huex_1_0_st_v2", values=rownames(frp4), mart=mart)

x<-read.table("/home/bianca/Downloads/up_combined_ids.txt")
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")

seq<-getSequence(id="CDH1", type="hgnc_symbol", seqType="cdna",mart=homo.anno)
seq<-getSequence(id=x[,1], type="ensembl_gene_id", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/bianca/Downloads/mrna_microara_bi_up.fasta")

micbio<-read.table("/home/bianca/Downloads/microarray_UP_prioritized.txt", sep="\t")
rnabio<-read.table("/home/bianca/Downloads/rnaseq_UP_prioritized.txt", sep="\t")
ridd_all<-read.table("/home/bianca/Downloads/mrna-microarray-out_ids_all", sep="\t")
ridd_exact<-read.table("/home/bianca/Downloads/mrna-microarray-exact-seq", sep="\t")


common_de_prio<-read.table("/home/bianca/Downloads/common_de_prio.tsv", sep="\t")
combined_de_prio<-read.table("/home/bianca/Downloads/Combined_de_prio.tsv", sep="\t")
# #annotation of genes
httr::set_config(httr::config(ssl_verifypeer = FALSE))



homo.anno = useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast")

list<-c("DH1", "CDK11A", "USP45","ZNF23")
seq<-getSequence(id=list, type="hgnc_symbol", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/alice_selected.fasta")

valuesr8<-rownames(rnaseq8)
valuesr24<-rownames(rnaseq24)

valuesm=micbio$V2
valuesr=rnabio$V2
values=c(micbio$V2, rnabio$V2)
values=c(common_de_prio$V2,combined_de_prio$V2)

attributes = c("ensembl_gene_id", "hgnc_symbol", "go_linkage_type")
attributes <- c("ensembl_gene_id", "hgnc_symbol","namespace_1003", "name_1006", "go_id", "go_linkage_type", "definition_1006")
attributes <- c("ensembl_gene_id", "hgnc_symbol","namespace_1003")
attributes <- c("ensembl_gene_id", "hgnc_symbol","name_1006", "go_id")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "definition_1006")
attributes <- c("ensembl_gene_id", "hgnc_symbol")


filters="ensembl_gene_id"
filters="hgnc_symbol"
values_ridd<-ridd_all$V1



# #listAttributes(homo.anno)
genes_ridd <- getBM(attributes = attributes, filters=filters,values = values_ridd, mart = homo.anno)


genes_ridd_prior <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

genescombined_common<-getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

genes_r8 <- getBM(attributes = attributes, filters=filters,values = valuesr8, mart = homo.anno)
genes_r24 <- getBM(attributes = attributes, filters=filters,values = valuesr24, mart = homo.anno)

#### prioritized combined and common genes that have ridd
now<-merge(genescombined_common,ridd_all, by.x=1, by.y=1)

###exact seq
now_Exact<-merge(genescombined_common,ridd_exact, by.x=1, by.y=1)

### bio info and exact and all ridd
########ensembl for vienna plots
grp<-merge(ridd_all, genes_ridd_prior, by.x=1, by.y=1)


ugrp<-unique(grp$hgnc_symbol)
grp1<-merge(ridd_exact, genes_ridd_prior, by.x=1, by.y=1)
ugrep1<-unique(grep1$ensembl_gene_id)


x1<-list(as.character(genes_ridd_prior$ensembl_gene_id),as.character(ridd_all$V1))
x2<-list(as.character(genes_ridd_prior$ensembl_gene_id),as.character(ridd_exact$V1))
x3<-list(as.character(micbio$V2),as.character(combined_de_prio$V2), as.character(micbio$V2),as.character(common_de_prio$V2))

a<-ggVennDiagram(
  x1, label_alpha = 0,
  category.names = c("Prioritized genes","Genes with CNGCNG")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")

b<-ggVennDiagram(
  x2, label_alpha = 0,
  category.names = c("Prioritized genes","Genes with CUGCAG")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")

c<-ggVennDiagram(
  x3, label_alpha = 0,
  category.names = c("Microarray","RNAseq","All UP regulated","Common up regulated")) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4") 

install.packages("gridExtra")       # Install gridExtra package
library("gridExtra") 


ridd_exact_m<-merge(ridd_exact, genesm, by.x = 1, by.y = 1)
urem<-unique(ridd_exact_m$hgnc_symbol)
ridd_all_m<-merge(ridd_all, genesm, by.x = 1, by.y = 1)
uram<-unique(ridd_all_m$hgnc_symbol)
dubr<-rnabio$V2[duplicated(rnabio$V2)]



genesr <- getBM(attributes = attributes, filters=filters,values = valuesr, mart = homo.anno)
ridd_exact_r<-merge(ridd_exact, genesr, by.x = 1, by.y = 1)
urer<-unique(ridd_exact_r$hgnc_symbol)
ridd_all_r<-merge(ridd_all, genesr, by.x = 1, by.y = 1)
urar<-unique(ridd_all_r$hgnc_symbol)
dubm<-micbio$V2[duplicated(micbio$V2)]


grp<-ridd_all[ridd_all$V1 %in% genes_ridd_prior$ensembl_gene_id,]
grp1<-merge(ridd_all, ridd_all, by.x=1, by.y=1)
