library(IsoformSwitchAnalyzeR)

sampleInfo<-data.frame(condition=c(rep(paste("D",8,sep=""),3),rep(paste("D",24,sep=""),3),rep(paste("M",8,sep=""),3),rep(paste("M",24,sep=""),3)),
                           sampleID=c(paste(rep(paste("D",8,sep=""),3),1:3,sep = "R"),paste(rep(paste("D",24,sep=""),3),1:3,sep="R"),paste(rep(paste("M",8,sep=""),3),1:3,sep="R"),paste(rep(paste("M",24,sep=""),3),1:3,sep="R")))
comp<-data.frame(condition_1=c("D24","D8"),
                 condition_2=c("M24","M8"))
stringtie<-importIsoformExpression("/home/bianca/Downloads/stringtie/finalbnew",readLength=100)

#myDesign<-model.matrix(~0 + condition,sampleInfo)

SwitchList <- importRdata(
      isoformCountMatrix   = stringtie$counts,
      isoformRepExpression = stringtie$abundance,
      designMatrix         = sampleInfo,
      isoformExonAnnoation = "/home/bianca/Downloads/stringtie/ens105.gtf",
      isoformNtFasta = "/home/bianca/Downloads/stringtie/ens105_transcriptome.fa",
      comparisonsToMake = comp,
      fixStringTieAnnotationProblem = TRUE,
      showProgress = FALSE,
      ignoreAfterBar = FALSE,ignoreAfterSpace = FALSE,
      ignoreAfterPeriod = TRUE
    )

backSwitch<-SwitchList
# ntSequence<-SwitchList$ntSequence
# aaSequence<-SwitchList$aaSequence
# exonsSequence<-SwitchList$exons

### Filter
SwitchList <- preFilter( SwitchList )

### Test for isoform switches
SwitchListDexseq <- isoformSwitchTestDEXSeq( SwitchList )   # OR
SwitchList <- isoformSwitchTestDRIMSeq( SwitchList ) 

### If analysing (some) novel isoforms (else use CDS from ORF as explained in importRdata() )
SwitchList <- addORFfromGTF( SwitchList )
SwitchList <- analyzeNovelIsoformORF( SwitchList )

### Extract Sequences
SwitchList <- extractSequence( SwitchList )


### Summary
extractSwitchSummary(SwitchList)    


#### confirmed isoforms
load("/home/bianca/Downloads/part1/isoform.RData")
genes<-as.data.frame(SwitchList$isoformFeatures$isoform_id)
genes$gene_name<-SwitchList$isoformFeatures$gene_name
colnames(genes)<-c("isoform","gene_name")

# homo.anno = useEnsembl(biomart = "ensembl", 
#                        dataset = "hsapiens_gene_ensembl", 
#                        mirror = "useast")
# 
# attributes <- c("ensembl_gene_id", "hgnc_symbol","ensembl_transcript_id")
# filters="ensembl_gene_id"
# 
# values<-SwitchList
# genes<- getBM(attributes = attributes, filters=filters,values = values, mart = homo.anno)


confirmed_exons_cor<-read.table("/home/bianca/Downloads/part1/all_confirmed_isoforms_ac",sep="\t")
exact<- confirmed_exons_cor[str_detect(confirmed_exons_cor$V3, "CUGCAG"),]

# genes$seq<-"NA"

# for (i in 1:nrow(genes)){
#   genes$seq[i]<-ifelse((genes[i,1] %in% confirmed_exons_cor$V1),"YES","NO")}
# 
# NO<-genes[genes$seq=="NO",]
# YES<-genes[genes$seq=="YES",]
# yes_and_no<-merge(YES,NO,by.x=2,by.y=2)
# ## POSSIBLY WRONG
# ridd_cor<-merge(yes_and_no,confirmed_exons_cor,by.x=2,by.y=1)

## merge confirmed exons with switchlist

is_it<- function(switch,con_exons){
  for (i in 1:nrow(switch)){
    switch$isoformFeatures$seq[i]<-ifelse((switch$isoformFeatures[i,3] %in% con_exons$V1),"YES","NO")
  }
  return(switch)
}

wow<-is_it(SwitchList, exact)

# de transcripts
de_trans<-SwitchList$isoformFeatures[(abs(SwitchList$isoformFeatures$dIF) > 0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05) ,]
# library(biomaRt)
# homo.anno = useEnsembl(biomart = "ensembl", 
#                        dataset = "hsapiens_gene_ensembl", 
#                        mirror = "useast")
# 
# attributes <- c("ensembl_gene_id", "hgnc_symbol","ensembl_transcript_id")
# filters="ensembl_transcript_id"
# seq<-getSequence(id=de_trans$isoform_id, type=filters, seqType="cdna",mart=homo.anno)
# 
# exportFASTA(seq, file="/home/bianca/Downloads/part1/elaborate_de_transcripts.fasta")

# switchlist that have fasta
#efasta<-de_trans[which(de_trans$isoform_id %in% seq$ensembl_transcript_id),]

#up 

up_down_seq<-function(SwitchList){
  a<-SwitchList$isoformFeatures[(SwitchList$isoformFeatures$dIF > 0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05) ,]
  upyes<-a[(a$seq=="YES"),]
  b<-SwitchList$isoformFeatures[(SwitchList$isoformFeatures$dIF < - 0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05) ,]
  updown<-b[(b$seq=="NO"),]
  omg<-merge(upyes,updown,by="gene_name")
  return(omg)
}

hey<-up_down_seq(wow)
# genes$hits<-"no_hit"
# for (i in 1:nrow(genes)){
#   genes$hits[i]<-ifelse((genes[i,1] %in% ridd_cor$ensembl_gene_id.x),"hit","no_hit")}
# 
# data <- data.frame(
#   group=c("hits","no hits"),
#   value=c(length(unique(genes[genes$hits=="hit",2])),length(unique(genes[genes$hits=="no_hit",2]))))
# data <- data %>% 
#   arrange(desc(group)) %>%
#   mutate(prop = value / sum(data$value) *100) %>%
#   mutate(ypos = cumsum(prop)- 0.5*prop )
# 
# data
# ggplot(data, aes(x="", y=value, fill=group)) +
#   geom_bar(stat="identity", width=1) +
#   coord_polar("y", start=0) +
#   theme_void() + 
#   theme(legend.position="right") +
#   #geom_text(aes(y = ypos, label = group), color = "black", size=6)  +
#   #scale_fill_brewer(palette="Set1") +
#   scale_colour_manual(values = c("skyblue4","palevioletred1")) +
#   #geom_text(aes(label = group), vjust = -1, nudge_y = 1, sie=6)  +
#   ggtitle("Up regulated 24 h") 
# 
# write.table(unique(genes[genes$hits=="hit",2]),'/home/ben/Desktop/exon_results/new/ids_common_hits_no_hits',sep="\t",quote=F,row.names=F)
# write.table(ridd_cor[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/24/isoforms_of_genes_24.txt",sep="\t",quote=F,row.names=F)
# write.table(unique(ridd_cor$hgnc_symbol),"/home/ben/Desktop/exon_results/new/24/for_bioinfominer_24.txt",sep="\t",quote=F,row.names=F)
# write.table(ridd_cor[,c(2,3,1,4,8,6,7)],"/home/ben/Desktop/exon_results/new/8/isoforms_of_genes_8.txt",sep="\t",quote=F,row.names=F)
# write.table(unique(ridd_cor$hgnc_symbol),"/home/ben/Desktop/exon_results/new/8/for_bioinfominer_8.txt",sep="\t",quote=F,row.names=F)
# 
# write.table(ridd_cor[,c(2,1,3,6,4,5)],"/home/bianca/Downloads/part1/ridd_possibly_all.txt",sep="\t",quote=F,row.names=F)
# 
# 



############### part 2
# test
# SwitchListAnalyzed <- analyzeCPAT(
#   switchAnalyzeRlist   = SwitchList,
#   pathToCPATresultFile = "/home/bianca/Downloads/part1/cpat.txt",
#   codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
#   removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
# )
# 
# 
# ### Add PFAM analysis
# SwitchListAnalyzed <- analyzePFAM(
#   switchAnalyzeRlist   = SwitchListAnalyzed,
#   pathToPFAMresultFile = "/home/bianca/Downloads/part1/pfam.txt",
#   showProgress=FALSE
# )
# 
# ### Add SignalP analysis
# SwitchListAnalyzed <- analyzeSignalP(
#   switchAnalyzeRlist       = SwitchListAnalyzed,
#   pathToSignalPresultFile  = "/home/bianca/Downloads/part1/signalIp.txt"
# )

### Add Iupred2a analysis
# SwitchListAnalyzed <- analyzeIUPred2A(
#   switchAnalyzeRlist        = SwitchListAnalyzed,
#   pathToIUPred2AresultFile = vec,
#   showProgress = FALSE
# )

# extractSwitchSummary(
#   SwitchListAnalyzed,
#   filterForConsequences = TRUE
# )
# hey<-extractTopSwitches(SwitchListAnalyzed, filterForConsequences = TRUE, n==Inf)
# 

# SwitchListAnalyzed
vec<-list.files("/home/bianca/Downloads/part1/website_results/", pattern="*.result", full.names  = TRUE)
SwitchList<-extractSequence(SwitchList)
# #all together
SwitchList <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = SwitchList,
  n                         = 10,
  removeNoncodinORFs        = TRUE,
  pathToCPATresultFile      =  "/home/bianca/Downloads/part1/cpat.txt",
  pathToPFAMresultFile      =  "/home/bianca/Downloads/part1/pfam.txt",
  pathToIUPred2AresultFile  = vec,
  pathToSignalPresultFile   = "/home/bianca/Downloads/part1/signalIp.txt",
  outputPlots               = FALSE,
  codingCutoff         = 0.725, 
)


extractTopSwitches(SwitchList, filterForConsequences = TRUE, n=Inf)

# ## if mega function was perfomrmed
# extractSwitchSummary(
#   SwitchList,
#   filterForConsequences = TRUE
# )

# check xbp1 for confirmation
options(scipen = 999)

# Down transcripts 
down_isoform_switch<-SwitchList$isoformFeatures[(SwitchList$isoformFeatures$dIF < -0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05),]
down_ridd_transcripts<-down_isoform_switch[up_isoform$gene_name %in% unique_gene_names,]
write.table(down_ridd_transcripts[,c(7,22,28)],"/home/bianca/Downloads/part1/bioinfo_DOWN_transcripts_confirmed.txt",sep="\t",quote = FALSE, row.names = FALSE)

up_isoform_switch<-SwitchList$isoformFeatures[(SwitchList$isoformFeatures$dIF > 0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05),]
up_ridd_transcripts<-up_isoform_switch[up_isoform$gene_name %in% unique_gene_names,]
up_ridd_genesNtrans<-up_ridd_transcripts[(up_ridd_transcripts$gene_log2_fold_change > 0.1 & up_ridd_transcripts$gene_switch_q_value < 0.05),]
write.table(up_ridd_transcripts[,c(7,22,28)],"/home/bianca/Downloads/part1/bioinfo_UP_transcripts_confirmed.txt",sep="\t",quote = FALSE, row.names = FALSE)
write.table(up_ridd_genesNtrans[,c(7,15,29)],"/home/bianca/Downloads/part1/bioinfo_UP_genes_up_transcripts_confirmed.txt",sep="\t",quote = FALSE, row.names = FALSE)


# UP GENES DE TRANSCRIPTS
de_isoform_switch<-SwitchList$isoformFeatures[(abs(SwitchList$isoformFeatures$dIF) > 0.1 & SwitchList$isoformFeatures$isoform_switch_q_value < 0.05),]
de_ridd_isoforms<-de_isoform_switch[de_isoform_switch$gene_name %in% unique_gene_names,]
up_ridd_genes_de_trans<-de_ridd_isoforms[(de_ridd_isoforms$gene_log2_fold_change > 0.1 & de_ridd_isoforms$gene_switch_q_value < 0.05),]
write.table(up_ridd_genes_de_trans[,c(7,15,29)],"/home/bianca/Downloads/part1/bioinfo_UP_genes_de_transcripts_confirmed.txt",sep="\t",quote = FALSE, row.names = FALSE)

bio<-read.table("/home/bianca/Downloads/part1/bioinfo_results_genes_up_trans_up.tsv", sep="\t")
bio<-read.table("/home/bianca/Downloads/part1/bioinfo_results_genes_DE_trans_up.tsv", sep="\t")
# wow <- function(i){
#   subset(
#     extractTopSwitches(SwitchList,filterForConsequences = TRUE,n=Inf,inEachComparison = TRUE)[,c('gene_name','condition_1','condition_2','gene_switch_q_value','Rank')],gene_name == i)
# }

unique_gene_names<-unique(ridd_cor$gene_name)
# 
# wowla<-function(unique_gene_names){
#   for (i in seq(1,length(unique_gene_names),1)){
#     wow(i)
#   }
# }

genes_that_exist<-unique_gene_names[unique_gene_names %in% SwitchList$switchConsequence$gene_name]

## unique genes from trans up and exact 
ugn<-unique(hey$gene_name)
wowla1<-function(unique_gene_names,SwitchList,c1,c2){
  for (i in seq(1,length(unique_gene_names),1)){
    switchPlot(SwitchList, gene = unique_gene_names[i],condition1 = c1,condition2 = c2)
  }
}
wowla1(ugn,wow,"D24","M24")
v2<-bio$V2
v2<-v2[-1]
wowla1(v2)
load("/home/bianca/Downloads/part1/isoform_part2_version2.RData")

## join switchlist with ridd cor

ridd_table<-ridd_cor %>% pivot_longer(cols=c(isoform.x,isoform.y),names_to = "trans",values_to="seq")


