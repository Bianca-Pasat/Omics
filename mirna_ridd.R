upmi1<-read.table("/home/ben/Documents/data/rankprod/on_cleaned/mirna/rp_mirna_mkc_vs_dmso_8_upregulatedBio.txt", sep="\t")
upmi2<-read.table("/home/ben/Documents/data/rankprod/on_cleaned/mirna/rp_mirna_mkc_vs_dmso_24_upregulatedBio.txt", sep="\t")
upmi<-rbind(upmi1,upmi2)
library(miRBaseConverter)
miRNANames=rownames(upmi)
version=checkMiRNAVersion(miRNANames,verbose = FALSE)
result=miRNA_NameToAccession(miRNANames,version=version)
Accessions = result$Accession
result1 = getMiRNASequence(Accessions)
result1<-na.omit(result1)

homo.anno = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")



filters="mirbase_id"
filters<-"mirbase_accession"
filters="ensembl_transcript_id"
attributes <- c("mirbase_id","mirbase_accession","ensembl_gene_id", "hgnc_symbol")
attributes <-c("ensembl_transcript_id","hgnc_id","entrezgene_id","refseq_mrna")
hgnc_trans_name
values<-c("ENST00000667629", "ENST00000366641","ENST00000449247","ENST00000442664")
values=rownames(upmi)
values=Accessions
genes_ridd <- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
genes_ridd <- getBM(attributes = attributes, filters=filters,values = rownames(upmi), mart = homo.anno)
filters="mirbase_accession"
attributes <- c("ensembl_gene_id", "hgnc_symbol")
genes_ridd <- getBM(attributes = attributes, filters=filters,values = Accessions, mart = homo.anno)

list<-c("DH1", "CDK11A", "USP45","ZNF23")
seq<-getSequence(id=list, type="hgnc_symbol", seqType="cdna",mart=homo.anno)


### yhc
yh<-read.table("/home/ben/Downloads/all_genes_mem_immune_related")
yh<-unique(yh)

m8<-read.table("/home/ben/Documents/data/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt", sep="\t")
m24<-read.table("/home/ben/Documents/data/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_24.txt", sep="\t")
values=rownames(m24)
filters="ensembl_gene_id"
attributes=c("ensembl_gene_id", "hgnc_symbol")
genes<- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
m24$ens<-rownames(m24)
m<-merge(m24,genes,by.x=3,by.y=1)
dim(m)

yh8<-m8[m8$V1 %in% yh$V1,]
yh8<-yh8[,c(1,2)]
colnames(yh8)<-c("names","FC_8")
yh24<-m[m$hgnc_symbol %in% yh$V1,]
yh24<-yh24[,c(4,2)]
colnames(yh24)<-c("names","FC_24")
yhf<-right_join(yh8,yh24,keep=TRUE)
write.table(yhf,"/home/ben/Desktop/yh_immune_cell_surface.txt", quote=F, row.names = F, sep="\t")