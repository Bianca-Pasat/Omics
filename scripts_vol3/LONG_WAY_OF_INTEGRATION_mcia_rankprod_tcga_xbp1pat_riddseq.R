demrna<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt", sep="\t")
demrna2<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_24.txt", sep="\t")

library(biomaRt)
ensembl<-useMart("ensembl")
ensembl<-useDataset("hsapiens_gene_ensembl", mart=ensembl)
filters="ensembl_gene_id"
attributes=c("ensembl_gene_id", "hgnc_symbol")
values=rownames(demrna2)
demrna2$names<-rownames(demrna2)
mrna_ids<-getBM(filters=filters, attributes=attributes, values=values, mart = ensembl)
demrna24<-merge(demrna2, mrna_ids, by.x="names", by.y="ensembl_gene_id")
demrna24_names<-demrna24[(demrna24$FC..class1.class2. > 0.2 | demrna24$FC..class1.class2. < - 0.2),]
demrna24_names01<-demrna24[(demrna24$FC..class1.class2. > 0.1 | demrna24$FC..class1.class2. < - 0.1),]
demrna802<-demrna[(demrna$V2 > 0.2 | demrna$V2 < -0.2), ]
demrna801<-demrna[(demrna$V2 > 0.1 | demrna$V2 < -0.1), ]

write.table(demrna24,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_bio_24_genenames", sep="\t", quote=F, row.names = T)
write.table(demrna24_names,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_bio_24_genenamesFC02", sep="\t", quote=F, row.names = T)
write.table(demrna24_names01,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_bio_24_genenamesFC01", sep="\t", quote=F, row.names = T)
write.table(demrna802,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_bio_8_genenamesFC02", sep="\t", quote=F, row.names = T)
write.table(demrna801,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_bio_8_genenamesFC01", sep="\t", quote=F, row.names = T)

dem01<-c(demrna24_names01$hgnc_symbol, rownames(demrna801))
dem02<-c(demrna24_names$hgnc_symbol, rownames(demrna802))
write.table(dem02,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_FC02", sep="\t", quote=F, row.names = T)
write.table(dem01,"/home/bianca/Downloads/integration_rankprod_mcia/de_mrna_FC01", sep="\t", quote=F, row.names = T)

demrnaA<-as.data.frame(unique.data.frame(demrna24[,c(4,2,3)]))
demrna24Nodup<-demrnaA[!duplicated(demrnaA$hgnc_symbol),]


demrna8Nodup<-demrna8[!duplicated(demrna8$hgnc_symbol),]
rownames(demrna8Nodup)<-demrna8Nodup$V1
demrna8Nodup<-demrna8Nodup[,-1]
colnames(demrna8Nodup)<-demrna8Nodup[1,]


##### proteins
deprot8<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/proteins/rp_clean_prot_mkc_vs_dmso_8_DEregulatedBio.txt", sep="\t")
deprot24<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/proteins/rp_clean_prot_mkc_vs_dmso_24_DEregulatedBio.txt", sep="\t")

deprot801<-deprot8[(deprot8$FC..class1.class2. > 0.1 | deprot8$FC..class1.class2. < -0.1 ),]
deprot2401<-deprot24[(deprot24$FC..class1.class2. > 0.1 | deprot24$FC..class1.class2. < -0.1 ),]
deprot802<-deprot8[(deprot8$FC..class1.class2. > 0.2 | deprot8$FC..class1.class2. < -0.2 ),]
deprot2402<-deprot24[(deprot24$FC..class1.class2. > 0.2 | deprot24$FC..class1.class2. < -0.2 ),]

write.table(deprot802,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_bio_8_genenamesFC02", sep="\t", quote=F, row.names = T)
write.table(deprot801,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_bio_8_genenamesFC01", sep="\t", quote=F, row.names = T)
write.table(deprot2402,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_bio_24_genenamesFC02", sep="\t", quote=F, row.names = T)
write.table(deprot2401,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_bio_24_genenamesFC01", sep="\t", quote=F, row.names = T)

dep01<-c(rownames(deprot2401), rownames(deprot801))
dep02<-c(rownames(deprot2402), rownames(deprot802))
write.table(dep02,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_FC02", sep="\t", quote=F, row.names = T)
write.table(dep01,"/home/bianca/Downloads/integration_rankprod_mcia/de_prot_FC01", sep="\t", quote=F, row.names = T)

################### mirnas
mirnas_names<-read.table("/media/bianca/My Passport/to_send_march/analysis_march23/de_rankprod_mirnas_ensembl.txt", sep="\t")
mirna8<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/mirna/rp_mirna_mkc_vs_dmso_8_DEregulatedBio.txt", sep="\t")
mirna24<-read.table("/home/bianca/Downloads/all_data_omics/OMICS_JANUARY_2023/rankprod/on_cleaned/mirna/rp_mirna_mkc_vs_dmso_24_DEregulatedBio.txt", sep="\t")

mirna801<-mirna8[(mirna8$FC..class1.class2. > 0.1 | mirna8$FC..class1.class2. < - 0.1),] ### same numbers
mirna802<-mirna8[(mirna8$FC..class1.class2. > 0.2 | mirna8$FC..class1.class2. < - 0.2),] ### only 6 out
mirna803<-mirna8[(mirna8$FC..class1.class2. > 0.3 | mirna8$FC..class1.class2. < - 0.3),]
mirna805<-mirna8[(mirna8$FC..class1.class2. > 0.5 | mirna8$FC..class1.class2. < - 0.5),]
mirna2401<-mirna24[(mirna24$FC..class1.class2. > 0.1 | mirna24$FC..class1.class2. < - 0.1),]
mirna2402<-mirna24[(mirna24$FC..class1.class2. > 0.2 | mirna24$FC..class1.class2. < - 0.2),]
mirna2403<-mirna24[(mirna24$FC..class1.class2. > 0.3 | mirna24$FC..class1.class2. < - 0.3),]
mirna2405<-mirna24[(mirna24$FC..class1.class2. > 0.5 | mirna24$FC..class1.class2. < - 0.5),]

demirnas0302<-c(rownames(mirna803), rownames(mirna2402))
demirnas0503<-c(rownames(mirna805), rownames(mirna2403))
demirnas05<-c(rownames(mirna805), rownames(mirna2405))
demirnas0302<-gsub("R", "r", demirnas0302)
demirnas0302<-gsub("-5p", "", demirnas0302)
demirnas0302<-gsub("-3p", "", demirnas0302)
demirnas0503<-gsub("R", "r", demirnas0503)
demirnas0503<-gsub("-5p", "", demirnas0503)
demirnas0503<-gsub("-3p", "", demirnas0503)
demirnas05<-gsub("R", "r", demirnas05)
demirnas05<-gsub("-5p", "", demirnas05)
demirnas05<-gsub("-3p", "", demirnas05)
demirnas0302_names<-mirnas_names[mirnas_names$V1 %in% demirnas0302,]
demirnas0503_names<-mirnas_names[mirnas_names$V1 %in% demirnas0503,]
demirnas05_names<-mirnas_names[mirnas_names$V1 %in% demirnas05,]


### lipids

########## intersect to see if you have less than 2000 for bio
mcia_integrated<-read.table("/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_allmirnas", sep="\t")


de01<-c(dep01, dem01)
de02<-c(dep02, dem02)

de020302<-c(dep02, dem02,demirnas0302_names$V3)
de020503<-c(dep02, dem02,demirnas0503_names$V3)
de0205<-c(dep02, dem02,demirnas05_names$V3)


mcia01<-mcia_integrated[mcia_integrated$V1 %in% de01,]
mcia02<-mcia_integrated[mcia_integrated$V1 %in% de02,]
mcia020302<-mcia_integrated[mcia_integrated$V1 %in% de020302,]
mcia020503<-mcia_integrated[mcia_integrated$V1 %in% de020503,]
mcia0205<-mcia_integrated[mcia_integrated$V1 %in% de0205,]

write.table(mcia01,"/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_01.txt", sep="\t", quote=F, row.names = F)
write.table(mcia02,"/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_02.txt", sep="\t", quote=F, row.names = F)
write.table(mcia020302,"/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_020302.txt", sep="\t", quote=F, row.names = F)
write.table(mcia020503,"/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_020503.txt", sep="\t", quote=F, row.names = F)
write.table(mcia0205,"/media/bianca/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_0205.txt", sep="\t", quote=F, row.names = F)


attributes = c("hgnc_symbol", "promoter_region_1000")
filters="hgnc_symbol"
values=mcia0205
promoters <- getBM(attributes = attributes, filters = filters, values = values, mart = ensembl)
# Define the length of the upstream and downstream regions
upstream_length <- 1000
downstream_length <- 0
attributes = c("hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position", "end_position", "transcription_start_site", "strand")
filters = "hgnc_symbol"
priori1<-read.table("/media/bianca/My Passport/to_send_march/gene_priori_mcia_fc05.tsv", sep="\t")
priori2<-read.table("/media/bianca/My Passport/to_send_march/gene_priori_first2000_integrated_nomirna.tsv", sep="\t")
values=c(priori1$V2, priori2$V2)

sign8<-demrna8[demrna8$V1 %in% values,]
sign24<-demrna24[demrna24$hgnc_symbol %in% values,]
signp8<-deprot8[rownames(deprot8) %in% values,]
signp24<-deprot24[rownames(deprot24)  %in% values,]


de_features_across_times<-full_join(sign8,sign24, signp8, signp24, keep=TRUE)

# Get the promoter sequences for the specified genes
promoters <- getBM(
  attributes = attributes,
  filters = filters,
  values = values,
  mart = ensembl
)

promoters1<-promoters[,c(1,3,4,5,7)]
promoters1<-unique.data.frame(promoters1)
promoters=promoters1
# Extract the promoter sequences from the reference genome
library(Biostrings)
promoter_sequences <- getSequence(
  id = values,
  type="hgnc_symbol",
  seqType="coding_gene_flank",
  upstream = 1000,
  mart=ensembl
)

library(stringr)
consensus_seq<-"GATGACGTG[TG][ATCG]{3}[AT]T"
consensus_seq<-"CCACG"
consensus_seq1<-"CGTCG"


for(i in 1:nrow(promoter_sequences)) {
  promoter_sequences$xbp1a[i] <- str_detect(promoter_sequences$coding_gene_flank[i], consensus_seq)
  promoter_sequences$xbp1b[i] <- str_detect(promoter_sequences$coding_gene_flank[i], consensus_seq1)
}
table(promoter_sequences$xbp1)
}


# promoter_sequences <- getSequence(
#   chromosome = promoters$chromosome_name,
#   start = promoters$start_position,
#   end = promoters$end_position,
#   type="hgnc_symbol",
#   seqType="cdna",
#   mart=ensembl
# )
# Add the promoter sequences to the data frame


# WE BLOCKED IRE1 SO XBP1 SIGN WILL BE DOWNREGULATED
xbp1_sign<-de_features[de_features$FC..class1.class2. < 0,]
ridd_sign<-de_features[de_features$FC..class1.class2. > 0,]
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")


GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)
#input SE,genesign,typesample and threshold
se <- data

# set your sample type e.g. typesample = "TP"
typesample <- c("TP", "NT") # (solid tumors)
typesample <- c("TP") # (solid tumors)
all<-colData(se)
tnbc<-all[which(all$paper_BRCA_Subtype_PAM50=="Basal"),]
se1<-se[,colnames(se) %in% tnbc$barcode]
exp <-se1[, TCGAquery_SampleTypes(colnames(se1), typesample=typesample)]

