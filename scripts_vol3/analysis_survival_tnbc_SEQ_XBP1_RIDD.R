library(Biostrings)
library(TCGAbiolinks)
library(biomaRt)
ensembl<-useMart("ensembl")
ensembl<-useDataset("hsapiens_gene_ensembl", mart=ensembl, mirror="uswest")
filters="ensembl_gene_id"
attributes=c("ensembl_gene_id", "hgnc_symbol")

mcia_integrated<-read.table("/media/ben/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_allmirnas.txt", sep="\t")
demrna8<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_8_GENE_NAMES.txt", sep="\t")
demrna2<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/mrna/rnaseq_rankprod_DEregulatedBio_24.txt", sep="\t")
library(biomaRt)
ensembl = useEnsembl(biomart = "ensembl", 
                       dataset = "hsapiens_gene_ensembl", 
                       mirror = "uswest")
ensembl = useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl", 
                     mirror = "useast")
filters="ensembl_gene_id"
attributes=c("ensembl_gene_id", "hgnc_symbol")
values=rownames(demrna2)
demrna2$names<-rownames(demrna2)
mrna_ids<-getBM(filters=filters, attributes=attributes, values=values, mart = ensembl)
demrna24<-merge(demrna2, mrna_ids, by.x="names", by.y="ensembl_gene_id")
colnames(demrna8)[1]<-"names"
colnames(demrna24)<-c("ens","FC..class1.class2.","P.value", "names")
test24<-demrna24[!(demrna24$names==""),]
test8<-demrna8[!demrna8$names=="",]

deprot8<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/proteins/rp_clean_prot_mkc_vs_dmso_8_DEregulatedBio.txt", sep="\t")
deprot24<-read.table("/home/ben/Desktop/OMICS_JANUARY_2023/rankprod/on_cleaned/proteins/rp_clean_prot_mkc_vs_dmso_24_DEregulatedBio.txt", sep="\t")
deprot8$names<-rownames(deprot8)
deprot24$names<-rownames(deprot24)
                         
all_demrnas<-full_join( test24,test8, by="names", keep=TRUE)
all_deprot<-full_join( deprot24,deprot8, by="names", keep=TRUE)

all_mrnas_mcia<-all_demrnas[(all_demrnas$names.x %in% mcia_integrated$V1 | all_demrnas$names.y %in% mcia_integrated$V1),]
all_deprot_mcia<-all_deprot[(all_deprot$names.x %in% mcia_integrated$V1 | all_deprot$names.y %in% mcia_integrated$V1),]
all_demrnas_mcia<-all_mrnas_mcia[,c(2,3,4,6,7,5)]
colnames(all_demrnas_mcia)=colnames(all_deprot_mcia)
all_FCvalues_mcia<-rbind(all_demrnas_mcia, all_deprot_mcia)

pp<-all_FCvalues_mcia
T1<-pp[(((pp$FC..class1.class2..x > 0 | pp$FC..class1.class2..x==NA ) & (pp$FC..class1.class2..y >0 | pp$FC..class1.class2..y==NA)) | ((pp$FC..class1.class2..x <0 | pp$FC..class1.class2..x==NA  ) & (pp$FC..class1.class2..y < 0 | pp$FC..class1.class2..y==NA ) )), ]

T1<-pp[(((pp$FC..class1.class2..x > 0 | pp$FC..class1.class2..x==NA ) & (pp$FC..class1.class2..y >0 | pp$FC..class1.class2..y==NA)) | ((pp$FC..class1.class2..x <0 | pp$FC..class1.class2..x==NA  ) & (pp$FC..class1.class2..y < 0 | pp$FC..class1.class2..y==NA ) )), ]

write.table(t1nona, "/media/ben/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_allmirnas_with_FCs", sep="\t", quote=F, row.names = T)
write.table(T1, "/media/ben/My Passport/to_send_march/MCIA-INTEGRATED/all_de_integrated_mcia_features_compund_names_allmirnas_with_FCs_AND_NANS", sep="\t", quote=F, row.names = T)

promoter_sequences<-list()
values=xbp1$names.x
promoter_sequence <- getSequence(
  id = values ,
  type="hgnc_symbol",
  seqType="coding_gene_flank",
  upstream = 1000,
  mart=ensembl
)

values=ridd$names.x
ridd_sequences <- getSequence(
  id = values,
  type="hgnc_symbol",
  seqType="cdna",
  mart=ensembl
)
exportFASTA(ridd_sequences, "/media/ben/My Passport/to_send_march/ridd_mcia.fasta")
exportFASTA(ridd_sequences, "/home/ben/Desktop/ridd_mcia.fasta")

exportFASTA(promoter_sequence, "/media/ben/My Passport/to_send_march/xbp1_mcia.fasta")
exportFASTA(promoter_sequence, "/home/ben/Desktop/xbp1_mcia.fasta")


############ survival analysis
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





#map and annotate gene symbols
#ensembl gene ids to external gene symbols
eg = as.data.frame(bitr(rownames(exp),
                        fromType="ENSEMBL",
                        toType="SYMBOL",
                        OrgDb="org.Hs.eg.db")) #drop=TRUE

eg <- eg[!duplicated(eg$SYMBOL),]#unique values
eg <- eg[!is.na(eg$SYMBOL),]#remove NA 


eg.ensembl <- as.character(eg$ENSEMBL)
dd <- match(eg.ensembl,rownames(exp))
dd <- na.omit(dd)
exp <- exp[dd, ]
rownames(exp) <- as.character(eg$SYMBOL)

#Print a message for genes of genesign were not found in exp.filt
not.found <- genesign[!(genesign %in% rownames(exp.filt))]
if(length(not.found) > 0){
  msg <- paste0("Sorry, I cant't find these genes: ", paste(not.found,collapse = " "))
  print(msg)
}

#keep the genes of signature that are present in exp.filt for further analysis
genesign <-  genesign[genesign %in% rownames(exp.filt)]
genesign.FC <- genesign.FC[rownames(genesign.FC) %in% genesign, ]


# preprocess



df_clean <- logCPM.genesign[rowSums(logCPM.genesign == 0) <= 10, ]
df_clean <- logCPM.genesign[apply(logCPM.genesign, 1, function(x) all(x >= 0)), ]

#create DGEList object (edgeR)
y.exp <- DGEList(counts = exp.filt)

#### Filtering for low counts BEFORE Normalization
#keep <- filterByExpr(y.exp)
#or 
keep <- rowSums(y.exp$counts==0) <= round(0.5*ncol(y))

y.exp <- y.exp[keep, , keep.lib.sizes=FALSE]

#TMM normalization
y.exp <- calcNormFactors(y.exp, method = "TMM")
#Note that by converting to CPMs we are normalising for the different sequencing depths for each sample.
cpm <- cpm(y.exp, log = FALSE, normalized.lib.sizes=TRUE)

# Filtering for low counts AFTER Normalization (optional...if you have filtered genes before normalization)
cpm.filt <- TCGAanalyze_Filtering(tabDF = cpm,
                                  method = "quantile",qnt.cut =  0.25)

# OR keep <- rowSums(cpm(y.exp)>0.5) >= round(0.5*ncol(y.exp))
#y.exp <- y.exp[keep, , keep.lib.sizes=FALSE]

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(cpm[,1],y.exp$counts[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#calculate log2CMP values
logCPM <- cpm(y.exp, prior.count=5, log=TRUE)

#map the genesign on logCPM.matrix
logCPM.genesign <- logCPM[genesign, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
q <- rowQuantiles(logCPM.genesign, probs=probs)
#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
rownames(aux) <- rownames(logCPM.genesign)
colnames(aux) <- colnames(logCPM.genesign)
#fill in the cells of aux matrix
for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(rownames(aux)[i] %in% XBP1sign) {   
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 1
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 2
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 3
      }else aux[i,j] <- 4
    }
    else if(rownames(aux)[i] %in% RIDDsign){
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 4
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 3
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 2
      }else aux[i,j] <- 1
    }
  }
}

#calculate a patient score
score <- colSums(aux)/length(genesign)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group <- "Mid expression"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group"] <- "Low expression"
colData(exp)[score >= high.cut,"level.group"] <- "High expression"


surv.data <- colData(exp)[colData(exp)$level.group %in% c("Low correlation", "High correlation"),]
clusterCol <- "level.group"

TCGAanalyze_survival(data = colData(exp),#data = surv.data,
                     clusterCol = clusterCol,
                     filename = NULL,#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)




### xbp1 and ridd separately
### XBP1

genesign<-XBP1sign
logCPM.genesign <- logCPM[genesign, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
q <- rowQuantiles(logCPM.genesign, probs=probs)
#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
rownames(aux) <- rownames(logCPM.genesign)
colnames(aux) <- colnames(logCPM.genesign)
#fill in the cells of aux matrix
for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(rownames(aux)[i] %in% XBP1sign) {   
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 1
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 2
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 3
      }else aux[i,j] <- 4
    }
  }
}

#calculate a patient score
score <- colSums(aux)/length(genesign)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group.xbp1 <- "Mid XBP1"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group.xbp1"] <- "Low XBP1"
colData(exp)[score >= high.cut,"level.group.xbp1"] <- "High XBP1"



##riddddd

genesign<-RIDDsign
logCPM.genesign <- logCPM[genesign, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# Row quantiles
q <- rowQuantiles(logCPM.genesign, probs=probs)
#create a scoring matrix assigning weights to the expression values of genes based on their expression level in the quantile distribution
aux <- matrix(rep(0), nrow=nrow(logCPM.genesign), ncol=ncol(logCPM.genesign))
rownames(aux) <- rownames(logCPM.genesign)
colnames(aux) <- colnames(logCPM.genesign)

for(j in 1:ncol(aux)){
  for(i in 1:nrow(aux)){
    if(rownames(aux)[i] %in% RIDDsign){
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 4
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 3
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 2
      }else aux[i,j] <- 1
    }
  }
}

#calculate a patient score
score <- colSums(aux)/length(genesign)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group.ridd <- "Mid RIDD"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group.ridd"] <- "Low RIDD"
colData(exp)[score >= high.cut,"level.group.ridd"] <- "High RIDD"


colData(exp)$grouping<-paste0(colData(exp)$level.group.xbp1, colData(exp)$level.group.ridd, sep="_")

surv.data <- colData(exp)[colData(exp)$level.group %in% c("Low correlation", "High correlation"),]
clusterCol <- "grouping"

TCGAanalyze_survival(data = colData(exp),#data = surv.data,
                     clusterCol = clusterCol,
                     filename = NULL,#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)

