###############This code is created by Konstantinos Voutetakis###############################################
#################################SurvivalPlot################################################################

# This module takes as input four arguments: a) a summarizedExperiment (SE) with TCGA gene expression and
# sample data with the following columns included:  days_to_death, vital_status, days_to_last_follow_up,   
# b) a gene signature (genesign) with gene_ids and fold changes, c) a survival threshold 
#############################################################################################################


library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

#input SE,genesign,typesample and threshold
load("H:/TCGA_GBM_TCGAbiolinks/TCGA_GBM.exp.rda")
se <- data
genesign <- read.csv("H:/TCGA_GBM_TCGAbiolinks/genesign.txt")
genesign <- as.character(genesign$genesign_ID)
genesign.FC <- as.matrix(genesign[,2])
rownames(genesign.FC) <- genesign$genesign_ID

# set your sample type e.g. typesample = "TP"
typesample <- "TP" # (solid tumors)


#Filtering SE for the user-defined sample type
exp <- se[, TCGAquery_SampleTypes(colnames(se), typesample=typesample)]

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

#create DGEList object (edgeR)
y.exp <- DGEList(counts = exp.filt)

#### Filtering for low counts BEFORE Normalization
#keep <- filterByExpr(y.exp)
#or keep <- rowSums(y.exp$counts==0) <= round(0.5*ncol(y))

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
    if(genesign.FC[i] > 0) {   
      if (logCPM.genesign[i,j] <= q[i,1]){ 
        aux[i,j] <- 1
      }else if (logCPM.genesign[i,j] > q[i,1] & logCPM.genesign[i,j] <= q[i,2]){
        aux[i,j] <- 2
      }else if(logCPM.genesign[i,j] > q[i,2] & logCPM.genesign[i,j] < q[i,3]){
        aux[i,j] <- 3
      }else aux[i,j] <- 4
    }
    else if(genesign.FC[i] < 0){
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

#survival Plot
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