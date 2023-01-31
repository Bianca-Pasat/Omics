library(ggVennDiagram)

mrna8 <- read.delim("//home/bianca/Desktop/RNAseq-avelogcpm_upper/8hours/rnaseq_rankprod_DEregulatedBio_8.txt", sep="", quote = "") 
mrna24 <- read.delim("/home/bianca/Desktop/RNAseq-avelogcpm_upper/24hours/rnaseq_rankprod_DEregulatedBio_24.txt", sep="", quote = "") 
mrna<-rbind(mrna8,mrna24)
mrna[,c(1,2)]<-as.numeric(mrna[,c(1,2)])

microarray8 <- read.delim("/home/bianca/Downloads/microarray_rankprod_DE_Bio_4.txt", sep="", quote = "") 
microarray24 <- read.delim("/home/bianca/Downloads/microarray_rankprod_DE_Bio_24.txt", sep="", quote = "") 
microarray<-rbind(microarray8,microarray24)


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
  ggplot2::scale_fill_gradient(low="skyblue3",high = "palevioletred2")

ggVennDiagram(
  x3, label_alpha = 0,
  category.names = c("RNAseq prioritized genes","Microarray prioritized genes")
) +
  ggplot2::scale_fill_gradient(low="deepskyblue4",high = "deeppink4")


upcom<-unique(data.frame(c(rownames(up_mrna), rownames(up_microarray))))

up_mrna2<-up_mrna[up_mrna$FC..class1.class2. > 0.15,]
up_mrna2<-up_mrna2[-c(1546:1570),]
up_microarray2<-up_microarray[up_microarray$FC..class1.class2. > 0.2,]
upcom2<-unique(data.frame(c(rownames(up_mrna2), rownames(up_microarray2))))

ridd_pred<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mrna-microarray-exact-seq", sep="\t")
ridd_pred2<-read.table("/home/ben/Desktop/miRNA_mRNA_thessaloniki_AND_microarray/mrna-microarray-out_ids_all", sep="\t")

ridd_found<-upcom[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred$V1,]
ridd_found_all<-ridd_pred2[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred2$V1,]

ridd_found_priori<-upcom[upcom$c.rownames.up_mrna...rownames.up_microarray.. %in% ridd_pred$V1,]

write.table(upcom, "/home/bianca/Downloads/up_all_combined_ids.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(upcom2, "/home/bianca/Downloads/up_all_combined_ids_02.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(rownames(common_up), "/home/bianca/Desktop/common_up_ids.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ridd_found, "/home/ben/Desktop/ridd_found.txt", sep="\t", quote=FALSE)

library(biomaRt)
homo.anno <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "hsapiens_gene_ensembl", 
                        mirror = "asia")
seq<-getSequence(id=x, type="ensembl_id", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/ben/Desktop/RIDD-EXON_USAGE_JANUARY_2022/second_analysis_exons_only/ridd/primary_inputs/ridd_input.fasta")
