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
library(sva)
library(tidyverse)
library("DEXSeq")
library("reshape2")
library("ggplot2")
library("stringr")
library(edgeR)
library(preprocessCore)
#setwd("/home/ben/Desktop/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/")
#path="/home/ben/Desktop/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq"

setwd("/home/bianca/Downloads/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/")
# path8="/home/bianca/Downloads/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/8"
# path24="/home/bianca/Downloads/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/24"
# pathgtf="/home/bianca/Downloads/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/"
path="/home/bianca/Downloads/rnaseq-aitor-feauturesCounts-DEXseq/DEXseq/"
countFiles = list.files(path, pattern=".bam$", full.names=TRUE)
#old
flattenedFile = list.files(path,pattern="gtf$", full.names=TRUE)

stringtie<-importIsoformExpression("/home/bianca/Downloads/stringtie/finalbnew",readLength=100)
# countFiles = list.files(path8, pattern=".bam$", full.names=TRUE)
# countFiles = list.files(pathgtf, pattern=".bam$", full.names=TRUE)
# 
# countFiles24 = list.files(path24, pattern=".bam$", full.names=TRUE)
### new with batch correction first
a<-read.table(countFiles[1],sep="\t")
countfiles<-countFiles[-1]
for (i in countfiles){
  la<-read.table(i, sep="\t")
  a<-merge(a, la, by=1)
}
# readin<- function(countFiles){
#   a<-read.table(countFiles[1],sep="\t")
#   countFiles<-countFiles[-1]
#   for (i in countFiles){
#     la<-read.table(i, sep="\t")
#     a<-merge(a, la, by=1)
#     return(a)
#     }
# }

#test
# hey<-read.table(countFiles[1],sep="\t")
# for (i in countfiles){
#   la<-read.table(i, sep="\t")
#   hey<-cbind(hey, la)
# }



sampleTable = data.frame(
  sample=c(paste("D24",1:3,sep = "R"),paste("D8",1:3,sep = "R"),paste("M24",1:3,sep = "R"),paste("M8",1:3,sep = "R")),
  treatment= c(rep("DMSO",6),rep("MKC",6)),
  libType = c(rep("single-end",12)),
  times=c( rep(c(rep(24,3), rep(8,3)),2)  ) ,
  replicate =as.factor(c(rep(c(1,2,3),4))))



counts<-a
counts<-a[-c(1:5),]
#colnames(counts)<-counts[length(counts[[1]]),]
#counts<-counts[-length(counts[[1]]),]
rownames(counts)<-counts$V1
counts=counts[,-1]
colnames(counts)<-sampleTable$sample
counts[,1:12]<-lapply(counts[1:12],as.numeric)

counts1<-counts

# keep.exprs <- aveLogCPM(counts) > 0
# counts <- counts[keep.exprs,]
# 
# ##upper quartile
# function.norm <- function(counts){
#   norm= normalize.quantiles.robust(as.matrix(counts), copy = TRUE) 
#   rownames(norm)=rownames(counts)
#   colnames(norm)=colnames(counts)
#   
#   return(norm)
# }


# norm_samples = function.norm(counts)
# 
# 
# 
# ### 8 hours
# norm_samples8 <-norm_samples[,c(4,5,6,10,11,12)]
# sampleTable8<-sampleTable[c(4,5,6,10,11,12),]

norm_samples8 <-counts1[,c(4,5,6,10,11,12)]
sampleTable8<-sampleTable[c(4,5,6,10,11,12),]

# 
# ###### 24 hours
# norm_samples24 <-norm_samples[,c(1,2,3,7,8,9)]
# sampleTable24<-sampleTable[c(1,2,3,7,8,9),]

norm_samples24 <-counts1[,c(1,2,3,7,8,9)]
sampleTable24<-sampleTable[c(1,2,3,7,8,9),]





#construct an DEXSeqDataSet object from this data
sampleData<-sampleTable8
colnames(sampleData)<-c("sample","condition","libtype","times","replicate")
cf<-countFiles[c(4,5,6,10,11,12)]
dxd1 = DEXSeqDataSetFromHTSeq(
  cf,
  sampleData=sampleData,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# exon name FIX THIS TO GET IT WITH STRSPLIT
fid<-dxd1@rowRanges$featureID
# gene name
gid<-dxd1@rowRanges$groupID

## rename this one to dxd
seqname<-dxd1@rowRanges@seqnames
fr<-dxd1@rowRanges@ranges
strand<-dxd1@rowRanges@strand
featureRa<-GRanges(seqname,fr,strand)
rm(dxd1)
### new dxd

sampleTable<-sampleTable24
sampleTable<-sampleTable8
sampleTable<-sampleTable24

poscleaned_count<-cleaned_count - min(cleaned_count)
counts<-round(poscleaned_count)
sampleData<-sampleTable
colnames(sampleData)<-c("sample","condition","libtype","times","replicate")


# 8 hours
counts=counts[,c(4,5,6,10,11,12)]
dxd=DEXSeqDataSet(
  counts,
  sampleData=sampleData,
  design= ~ sample + exon + condition:exon,
  featureID= fid,
  groupID = gid,
  featureRanges = featureRa)


## with cfs
cf<-countFiles[c(1,2,3,7,8,9)]
dxd = DEXSeqDataSetFromHTSeq(
  cf,
  sampleData=sampleData,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
# first 6 columns = number of reads mapping to exonic regions, last 6 = sum of the counts for rest of the exons from the same gene on each sample.

### Normalization (same method as in DEseq2)
dxd = estimateSizeFactors( dxd )

# 8 
dxd1 = estimateSizeFactors( dxd1 )

# Dispersion estimation (estimate strength of the noise)
dxd = estimateDispersions( dxd )
dxd1 = estimateDispersions( dxd1 )
#plotDispEsts( dxd )

### Testing for differential exon usage
# fits a generalized linear model with the formula ~sample + exon + condition:exon and compare it to the smaller model (the null model) ~ sample + exon.

dxd = testForDEU( dxd )
dxd1 = testForDEU( dxd1 )

# estimate relative exon usage fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

dxd1 = estimateExonFoldChanges( dxd1, fitExpToVar="condition")
### Results summary
dxr1 = DEXSeqResults( dxd )
mcols(dxr1)$description
resultados <- as.data.frame(dxr1)

# 8 with cfs no sva
dxr2 = DEXSeqResults( dxd1 )
mcols(dxr2)$description
resultados2 <- as.data.frame(dxr2)

#RENAME TREATMENT TO CONDIITION
### VISUALIZATION
# MA plot
#plotMA( dxr1, cex=0.8 )
#table ( dxr1$padj < 0.05 )
hits <- subset(resultados, resultados$padj < 0.05)
hitsire1<-subset(hits, hits$log2fold_MKC_DMSO > 0 )
hitsire1WOW<-subset(hits, hits$log2fold_MKC_DMSO > 1 )

hits_table <- as.data.frame(table(hitsire1$groupID))
hits_tablefeauture <- as.data.frame(table(hitsire1$groupID))
#hits_table <- as.data.frame(table(hitsire1$groupID))

dnew<-subset(dxr1, dxr1$padj < 0.05)

# plot exon usage
dexplot<-function(ddx){
  for (i in seq(1,length(hitsire1$groupID),1)){
    a<-hitsire1$groupID[i]
    plotDEXSeq(ddx, a, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
  }
}

plotDEXSeq(dxr1, "ENSG00000141736", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1, "ENSG00000148396", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


# load results from previous analysis
load("/home/bianca/Desktop/DEXSEQ_exons_results/exons_results/new/exon_analysis_8.RData")
plotDEXSeq(dxr1,"ENSG00000163681", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000214022", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000176014", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
load("/home/bianca/Desktop/DEXSEQ_exons_results/exons_results/new/exon_analysis_24.RData")
plotDEXSeq(dxr1,"ENSG00000163681", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000214022", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000176014", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# 8 hours with cfs SLMAP, REPIN1, TUBB6 
plotDEXSeq(dxr2,"ENSG00000163681", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr2,"ENSG00000214022", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr2,"ENSG00000176014", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

# 24 hours with cfs SLMAP, REPIN1, TUBB6 
plotDEXSeq(dxr1,"ENSG00000163681", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000214022", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq(dxr1,"ENSG00000176014", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#dexplot(dnew)
#dexplot(dxr1)

### Save summary results in HTML bundle
#DEXSeqHTML(dnew, FDR = 0.01)

#setwd("/home/bianca/Desktop/exon_results/24")
setwd("/home/bianca/Desktop/exon_results/8")
saveRDS(hitsire1,"/home/bianca/Desktop/hitsire124hours")

write.table(hitsire1,"/home/bianca/Desktop/hitsire18hours.txt",quote=F, row.names = F)

DEXSeqHTML(dxr1, FDR = 0.01)

#plot frecuencies pie
exon_frequency_pie <- function(nombre_gen) {
  nombre_gen1  <- subset(resultados, resultados$groupID == nombre_gen)
  m1_sums <- na.omit(nombre_gen1$DMSO)/log2(na.omit(nombre_gen1$genomicData.width))
  m2_sums <- na.omit(nombre_gen1$MKC)/log2(na.omit(nombre_gen1$genomicData.width))
  par(mfrow = c(1,2))
  pie(m1_sums, main="DMSO")
  pie(m2_sums, main="MKC")
}  

exon_frequency_pie("ENSG00000048162")

# ### GET SEQUENCES OF THE EXONS
# ### add another column with the nucleotide seqs
# ## might want to look at range around it, not just there
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# 
# 
# 
# chromO <- hitsire1$genomicData.seqnames
# starT <-  hitsire1$genomicData.start
# enD <- hitsire1$genomicData.end
# 
# 
# sequence_exon <- function(chrom0, starT, enD){
#   for (i in seq(1, length(chrom0))){
#     b<-as.character(getSeq(Hsapiens, chrom0[i], start = starT[i] - 100, end = enD[i] + 100))
#   }
# }
# 
# fa<-sequence_exon(chromO,starT,enD)
# 
# exportFASTA(fa, file="/home/bianca/Desktop/my_exons24.fasta")
# exportFASTA(fa, file="/home/bianca/Desktop/my_exons8.fasta")
# 
# 
# #aitor
# hits_table_all_annot <- read.csv("C:/Users/alman/Desktop/Grecia/Exon_usage/hits_table_all_annot.csv")
# 
# sequence_exon <- function(x){
#   as.character(getSeq(Hsapiens, chromO[x], start = starT[x], end = enD[x]))
# }
# 
# sequence_exon_flank <- function(x) {
#   as.character(getSeq(Hsapiens, chromO[x], start = starT[x]-100, end = enD[x]+100))
# }
# 
# for(i in 1:nrow(hits)) {
#   hits$sequence[i] <- sequence_exon(i)
# }
# 
# for(i in 1:nrow(hits)) {
#   hits$sequence_flanks[i] <- sequence_exon_flank(i)
# }
# 
# for(i in 1:nrow(hitsire1)) {
#   hitsire1$sequence_flanks[i] <- sequence_exon_flank(i)
# }
# 
# ### scan for consensus cleavage sequence in exon seqs
# consensus_seq <- as.character("CTGCAG")
# 
# for(i in 1:nrow(hits)) {
#   hits$cleavage[i] <- str_detect(hits$sequence[i], consensus_seq)
# }
# 
# 
# 
# for(i in 1:nrow(hitsire1)) {
#   hitsire1$cleavage_flank[i] <- str_detect(hitsire1$sequence_flanks[i], consensus_seq)
# }
# 
# #### make fasta from columns
# ire1_exon_seq<-hitsire1[which(hitsire1$cleavage_flank == "TRUE"),]
# ire1_exon_seq1<-ire1_exon_seq[,c("groupID", "sequence_flanks")]
# write.table(ire1_exon_seq1,"/home/bianca/Desktop/exon_results/8/my_exons_seq.txt", sep="\t")
# write.table(ire1_exon_seq1,"/home/bianca/Desktop/exon_results/24/my_exons_seq.txt", sep="\t")
# 
# ############################################
# ### Screen ALL exons for a given gene
# ### TAKE A LOOK AT EACH GENE AT A TIME
# library("biomaRt")

# listEnsembl(GRCh=37)
# listEnsembl(version=99)
# 
# ensembl = useMart("ensembl")
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# 
# # LOOK FOR CONSENSUS CLEAVAGE SEQUENCE trough all exons of one gene
# CD44_exons <- getSequence( id="CD44" ,type='hgnc_symbol',seqType = 'gene_exon', mart = ensembl)
# # CD44 <- getSequence(chromosome = 4, start = 50000-3000, end = 50000+3000,type='hgnc_symbol',seqType = 'gene_exon_intron', mart = GENES)
# 
# for(i in 1:nrow(CD44_exons)) {
#   CD44_exons$cleavage[i] <- str_detect(CD44_exons$gene_exon[i], consensus_seq)
# }
# 
# 
# ### ANY EXON WITH CLEAVAGE IN GENE X?
# cleavage_true <- function(x) {
#   exons <- getSequence( id= x ,type='hgnc_symbol',seqType = 'gene_exon', mart = ensembl)
#   for(i in 1:nrow(exons)) {
#     exons$cleavage[i] <- str_detect(exons$gene_exon[i], consensus_seq)
#   }
#   table(exons$cleavage)
# }
# 
# cleavage_true("CD44")
# 
# #FALSE  TRUE 
# #75    22 
# 
# for(i in 1:nrow(hits_table_all_annot)) {
#   hits_table_all_annot$cleavage[i] <- cleavage_true(hits_table_all_annot$hgnc_symbol[i])[2]
# }
# 
# ###
# CD44_exons <- subset(resultados, resultados$groupID == "ENSG00000026508")
# 
# for(i in 1:nrow(CD44_exons)) {
#   CD44_exons$sequence[i] <- sequence_exon(i)
# }
# 
# for(i in 1:nrow(CD44_exons)) {
#   CD44_exons$cleavage[i] <- str_detect(CD44_exons$sequence[i], consensus_seq)
# }
# 
# ### annotate gene symbol
# 
# probando <- function(x) {
#   getSequence(chromosome = hits_seq_cleav_annot$genomicData.seqnames[x], start = hits_seq_cleav_annot$genomicData.start[i], end = hits_seq_cleav_annot$genomicData.end[i],type='hgnc_symbol',seqType = 'gene_exon_intron', mart = ensembl)[2]
# }
# 
# for(i in 1:nrow(hits_seq_cleav_annot)) {
#   hits_seq_cleav_annot$hgnc_symbol[i] <- probando(i)
# }


homo.anno = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "useast")
ire1_exon_seq1<-c(hitsire1$groupID,hitsire1$featureID)
attributes <- c("ensembl_gene_id", "hgnc_symbol","ensembl_transcript_id")
filters="ensembl_gene_id"
two<-as.data.frame(do.call(rbind,str_split(hitsire1$groupID, "\\+")))
three<-rbind(two$V1,two$V2,two$V3,two$V4,two$V5,two$V6,two$V7)
values<-two$V1
genes<- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)
seq<-getSequence(id=genes$ensembl_transcript_id, type="ensembl_transcript_id", seqType="cdna",mart=homo.anno)
exportFASTA(seq, file="/home/bianca/Desktop/exons_results/24/exons_transcripts_24.fasta")



confirmed_exons_cor<-read.table("/home/ben/Desktop/exon_results/new/8/out_column_8hours_cor",sep=" ")
confirmed_exons_cor<-read.table("/home/ben/Desktop/exon_results/new/24/out_column_24hours_cor",sep=" ")

confirmed_exons_cor8<-read.table("/home/ben/Desktop/exon_results/new/8/out_column_8hours_cor",sep=" ")
confirmed_exons_cor24<-read.table("/home/ben/Desktop/exon_results/new/24/out_column_24hours_cor",sep=" ")
confirmed_exons_cor<-rbind(confirmed_exons_cor8,confirmed_exons_cor24)

confirmed_exons8<-read.table("/home/ben/Desktop/exon_results//new/8/out_column_8hours",sep="\t")
confirmed_exons24<-read.table("/home/ben/Desktop/exon_results/new/24/out_column_24hours",sep="\t")

genes$seq<-"NA"
for (i in 1:nrow(genes)){
  genes$seq[i]<-ifelse((genes[i,3] %in% confirmed_exons$V1),"YES","NO")}
for (i in 1:nrow(genes)){
  genes$seq[i]<-ifelse((genes[i,3] %in% confirmed_exons_cor$V1),"YES","NO")}
for (i in 1:nrow(genes)){
  genes$seq[i]<-ifelse((genes[i,3] %in% confirmed_exons_cor8$V1),"YES","NO")}
for (i in 1:nrow(genes)){
  genes$seq[i]<-ifelse((genes[i,3] %in% confirmed_exons_cor24$V1),"YES","NO")}
NO<-genes[genes$seq=="NO",]
YES<-genes[genes$seq=="YES",]
yes_and_no<-merge(YES,NO,by.x=2,by.y=2)
ridd_cor<-merge(yes_and_no,confirmed_exons_cor,by.x=3,by.y=1)
ridd_cor<-merge(yes_and_no,confirmed_exons_cor24,by.x=3,by.y=1)
ridd_cor<-merge(yes_and_no,confirmed_exons_cor8,by.x=3,by.y=1)
ridd_cor<-merge(yes_and_no,confirmed_exons,by.x=3,by.y=1)

genes$hits<-"no_hit"
for (i in 1:nrow(genes)){
  genes$hits[i]<-ifelse((genes[i,1] %in% ridd_cor$ensembl_gene_id.x),"hit","no_hit")}

data <- data.frame(
  group=c("hits","no hits"),
  value=c(length(unique(genes[genes$hits=="hit",2])),length(unique(genes[genes$hits=="no_hit",2]))))
data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

data
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  #geom_text(aes(y = ypos, label = group), color = "black", size=6)  +
  #scale_fill_brewer(palette="Set1") +
  scale_colour_manual(values = c("skyblue4","palevioletred1")) +
  #geom_text(aes(label = group), vjust = -1, nudge_y = 1, sie=6)  +
  ggtitle("Up regulated 24 h") 

write.table(unique(genes[genes$hits=="hit",2]),'/home/ben/Desktop/exon_results/new/ids_common_hits_no_hits',sep="\t",quote=F,row.names=F)
write.table(ridd_cor[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/24/isoforms_of_genes_24.txt",sep="\t",quote=F,row.names=F)
write.table(unique(ridd_cor$hgnc_symbol),"/home/ben/Desktop/exon_results/new/24/for_bioinfominer_24.txt",sep="\t",quote=F,row.names=F)
write.table(ridd_cor[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/8/isoforms_of_genes_8.txt",sep="\t",quote=F,row.names=F)
write.table(unique(ridd_cor$hgnc_symbol),"/home/ben/Desktop/exon_results/new/8/for_bioinfominer_8.txt",sep="\t",quote=F,row.names=F)

write.table(ridd_cor[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/common_isoforms_of_genes.txt",sep="\t",quote=F,row.names=F)


prioritized_8<-read.table("/home/ben/Desktop/exon_results/new/8/bioinfo",sep="\t")
prioritized_24<-read.table("/home/ben/Desktop/exon_results/new/24/bioinfo",sep="\t")

prio_iso<-ridd_cor[ridd_cor$hgnc_symbol %in% prioritized_8$V2,]
prio_iso<-ridd_cor[ridd_cor$hgnc_symbol %in% prioritized_24$V2,]

write.table(prio_iso[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/8/prioriizd_isoforms_of_genes_8.txt",sep="\t",quote=F,row.names=F)
write.table(prio_iso[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/24/prioriizd_isoforms_of_genes_24.txt",sep="\t",quote=F,row.names=F)


###
# library("biomaRt")
# human = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org", path="/biomart/martservice",ensemblRedirect = FALSE)
# 
# transcript_info = data.frame(unique(getBM(attributes = c("chromosome_name", "genomic_coding_start","genomic_coding_end"),filters="ensembl_transcript_id", values="ENST00000269305",mart = human)))
# 
# transcript_info[order(transcript_info$genomic_coding_start),]
# 
# ensembl_gene_id