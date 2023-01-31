library(ggVennDiagram)

mrna8 <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/mRNA/DE/8HOURS/CORRECT-ave_log_cpm-sva-combat/rnaseq_rankprod_DEregulatedBio_8.txt", sep="", quote = "") 
mrna24 <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/mRNA/DE/24HOURS/avelogcpm-sva-combat/rnaseq_rankprod_DEregulatedBio_24.txt", sep="", quote = "") 
mrna<-rbind(mrna8,mrna24)
mrna[,c(1,2)]<-as.numeric(mrna[,c(1,2)])

microarray8 <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/microarray_exp/MY_ANALYSIS-microarray_results_rankprod_sva_combat/4_hours_MKC_vs_DMSO/microarray_rankprod_DE_Bio_4.txt", sep="", quote = "") 
microarray24 <- read.delim("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/microarray_exp/MY_ANALYSIS-microarray_results_rankprod_sva_combat/24_hours_MKC_vs_DMSO/microarray_rankprod_DE_Bio_24.txt", sep="", quote = "") 
microarray<-rbind(microarray8,microarray24)

#proteomics-rnaseq
deprot8<-read.table("/home/bianca/Desktop/proteomics_de_results/rp_prot_mkc_vs_dmso_8_DEregulatedBio.txt",sep="\t")
deprot24<-read.table("/home/bianca/Desktop/proteomics_de_results/rp_prot_mkc_vs_dmso_24_DEregulatedBio.txt",sep="\t")

derna8<-read.table("/home/bianca/Downloads/ANALYSIS_APRIL/DE/RNAseq/rnaseq_rankprod_DEregulatedBio_8.txt")
derna24<-read.table("/home/bianca/Downloads/ANALYSIS_APRIL/DE/RNAseq/rnaseq_rankprod_DEregulatedBio_24.txt")

httr::set_config(httr::config(ssl_verifypeer = FALSE))
homo.anno = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")

values<-rownames(derna8)
values<-rownames(derna24)
attributes <- c("ensembl_gene_id", "hgnc_symbol")
filters="ensembl_gene_id"
dernasymbol8 <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
dernasymbol24 <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)

x<-list(rownames(deprot8),unique(dernasymbol8$hgnc_symbol))
x<-list(rownames(deprot24),dernasymbol24$hgnc_symbol)

com8<-dernasymbol8[dernasymbol8$hgnc_symbol %in% rownames(deprot8),]
com24<-dernasymbol8[dernasymbol24$hgnc_symbol %in% rownames(deprot24),]

write.table(com8$hgnc_symbol,"/home/bianca/Desktop/proteomics_de_results/com_ven8.txt",sep="\t",quote=F, row.names = F)
write.table(com24$hgnc_symbol,"/home/bianca/Desktop/proteomics_de_results/com_ven24.txt",sep="\t",quote=F, row.names = F)


ggVennDiagram(
  x, label_alpha = 0,
  category.names = c("Proteomics",'RNAseq')
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")


#########
up_mrna<-mrna[mrna[,1] > 0,]
down_mrna<-mrna[mrna[,1] < 0,]

up_microarray<-microarray[microarray[,1] > 0,]
down_microarray<-microarray[microarray[,1] < 0,]

bio_up_mrna<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mRNA-miRNA-thes/mRNA/DE/rnaseq_UP_prioritized.txt", sep="\t")[2]
bio_up_mrna<-bio_up_mrna[-1,]
bio_up_microarray<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/microarray_exp/MY_ANALYSIS-microarray_results_rankprod_sva_combat/microarray_UP_prioritized.txt", sep="\t")[2]
bio_up_microarray<-bio_up_microarray[-1,]
bio_common<-data.frame(c(bio_up_mrna,bio_up_microarray))

x<-list(rownames(mrna),rownames(microarray))
common_de<-mrna[rownames(mrna) %in% rownames(microarray),]

x1<-list(rownames(up_mrna),rownames(up_microarray))
common_up<-up_mrna[rownames(up_mrna) %in% rownames(up_microarray),]

x3<-list(as.character(bio_up_mrna),as.character(bio_up_microarray))



ggVennDiagram(
  x, label_alpha = 0,
  category.names = c("RNAseq","Microarray")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")



a<-ggVennDiagram(
  x1, label_alpha = 0,
  category.names = c("up RNAseq","up Microarray")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")

ggVennDiagram(
  x3, label_alpha = 0,
  category.names = c("RNAseq prioritized genes","Microarray prioritized genes")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")


upcom<-unique(data.frame(c(rownames(up_mrna), rownames(up_microarray))))
ridd_pred<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mrna-microarray-exact-seq", sep="\t")
ridd_pred2<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mrna-microarray-out_ids_all", sep="\t")

ridd_found<-upcom[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred$V1,]
ridd_found_all<-ridd_pred2[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred2$V1,]

ridd_found_priori<-upcom[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred$V1,]

write.table(upcom, "/home/ben/Desktop/up_combined_ids.txt", sep="\t", quote=FALSE)
write.table(rownames(common_up), "/home/ben/Desktop/common_up_ids.txt", sep="\t", quote=FALSE)
write.table(ridd_found, "/home/ben/Desktop/ridd_found.txt", sep="\t", quote=FALSE)

library(biomaRt)
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")
seq<-getSequence(id=x, type="ensembl_id", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/ridd/primary_inputs/ridd_input.fasta")
