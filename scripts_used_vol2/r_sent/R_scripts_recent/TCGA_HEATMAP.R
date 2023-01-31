#Rscript TCGA_HEATMAP.R -c 'TCGA-COAD' -g '/Users/user/NUIG-HCT116/wptb_comon_network/Gene_Prioritization.tsv' -R '/Users/user/hout.html' -o '/Users/user/hout.pdf'
rm(list = ls()) 

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler)) #
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(circlize)) 
suppressPackageStartupMessages(library(ComplexHeatmap)) #
suppressPackageStartupMessages(library(RankProd))
suppressPackageStartupMessages(library(getopt))

args = commandArgs(trailingOnly=TRUE)
BiocManager::install(c("ggraph","enrichplot","clusterProfiler"))

#cols=eval(parse(text=args[1]))
# cancer<-as.character(args[1])
# genesign<-read.table(args[2], header = T,  sep="\t")
# htmlPath<-args[3]
# outPath<-args[4]

# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "cancer", "c", 1, "character",
  "genesign", "g", 1, "character",
  "htmlPath", "R", 2, "character",
  "outPath", "o", 2, "character"),
byrow=TRUE, ncol=4)
opt <- getopt(spec)




# cancer<-as.character(args[1])
# genesign<-read.table(args[2], header = T,  sep="\t")
# cancer<-as.character(args[1])
# genesign<-read.table(args[2], header = T,  sep="\t")
# 

# Generate output folder and paths
makeOut <- function(filename) {
  return(paste0(opt$outPath, "/", filename))
}


###these can be put in matrix, if you have many comparisons 
heatOutPdf <- character()
heatOutPdf <- makeOut(paste0("heatmap", ".pdf"))



# Initialise data for html links and images, data frame with columns Label and
# Link
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
                        stringsAsFactors=FALSE)


	
 	
 	
 	
 	

###################################################################################################
##### HEATMAP #######################################
####################################################################################################

########### Download TCGA DATA
query <- GDCquery(project = opt$cancer, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")

query <- GDCquery(project = "TCGA-UVM", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")


GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)


# load("/Users/user/GDC_data/TCGA-ACC")
# GDCdownload(query, method = "client", files.per.chunk = 10)
# data <- GDCprepare(query)

se <- data

# set your sample type e.g. typesample = "TP"
# typesample <- c("TP", "NT") # (all tumors for difexp)
typesample <- c("TP") # (solid tumors)
exp <-se[, TCGAquery_SampleTypes(colnames(se), typesample=typesample)]
#Filtering SE for the user-defined sample type

# expTP <- TCGAquery_SampleTypes(barcode=colnames(exp), typesample = "TP")
# 
# expNT <- TCGAquery_SampleTypes(barcode=colnames(exp),typesample = "NT")

eg = as.data.frame(bitr(rownames(exp),
                        fromType="ENSEMBL",
                        toType="SYMBOL",
                        OrgDb="org.Hs.eg.db")) #drop=TRUE

eg <- eg[!duplicated(eg$SYMBOL),]#unique values
eg <- eg[!is.na(eg$SYMBOL),]#remove NA


eg.ensembl <- as.character(eg$ENSEMBL)
dd <- match(eg.ensembl,rownames(exp))

exp <- exp[dd, ]

rownames(exp) <- as.character(eg$SYMBOL)

# Filtering for low counts
exp.filt <- TCGAanalyze_Filtering(tabDF = assay(exp),
                                  method = "quantile",qnt.cut =  0.25)

#create DGEList object (edgeR)
y.exp <- DGEList(counts = exp.filt)

#TMM normalization
y.exp <- calcNormFactors(y.exp, method = "TMM")
###LOGCPM
logCPM <- cpm(y.exp, prior.count=5, log=TRUE)
clinical_patient_Cancer<-colData(se)
save.image("/Users/user/GDCdata/TCGA-UVM.RData")
#input SE,genesign,typesample and threshold
# load("H:/TCGA_GBM_TCGAbiolinks/TCGA_GBM.exp.rda")
# se <- data
#genesign <- read.csv("H:/TCGA_GBM_TCGAbiolinks/genesign.txt")
#genesign<-read.table("/Users/user/NUIG-HCT116/wptb_comon_network/Gene_Prioritization.tsv",  header = T,  sep="\t")
genesign<-read.table(opt$genesign, header=T,sep="\t")

genesign <- as.character(genesign[,2])

genesign.FC <- as.matrix(genesign)
rownames(genesign.FC) <- genesign


#keep the genes of signature that are present in exp.filt for further analysis
genesign <-  genesign[genesign %in% rownames(exp.filt)]
genesign.FC <- genesign.FC[rownames(genesign.FC) %in% genesign, ]

#map the genesign on logCPM.matrix
logCPM.genesign <- logCPM[genesign.FC, ]

zscores<-scale(logCPM.genesign)
test<-zscores



col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5

clinical_patient_Cancer<-clinical_patient_Cancer[rownames(clinical_patient_Cancer) %in% colnames(test),]

signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(tissue_type=clinical_patient_Cancer[,4])

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
pdf(heatOutPdf)
Heatmap(test, name="expression",
            col=col_fun,
            cluster_rows=T,
            cluster_columns = T,
            row_names_side = "left",
            show_row_names = T,
            #column_order = order(clinical_patient_Cancer$shortLetterCode),
            #row_order = order(rownames(logCPM.genesign)),
            show_column_names = F,
            #column_names_gp=gpar(cex=fontsize),
            show_row_dend = FALSE,
            show_column_dend=TRUE,
            row_names_gp=gpar(cex=fontsize+0.1),
            row_names_max_width = unit(5, "cm"),
            clustering_distance_rows ="euclidean",
            clustering_method_rows = "ward.D",
            clustering_distance_columns =  "euclidean",
            clustering_method_columns = "ward.D",
            row_dend_width = unit(10, "mm"),
            #row_km = 2,
            #column_km =2,
            #km=2,
            left_annotation=row,
            bottom_annotation = col,
            heatmap_width = unit(20, "cm"),
            width = NULL,
            heatmap_height = unit(20, "cm"))
draw(hm)
linkName <- paste0("heatmap", ".pdf")
linkAddr <- paste0("heatmap", ".pdf")
linkData <- rbind(linkData, c(linkName, linkAddr))
imageData<- rbind(linkData, c(linkName, linkAddr))
invisible(dev.off())
################################################################################
### HTML Generation
################################################################################
# Create cata function: default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)

cat("", file=opt$htmlPath)
cata <- function(..., file = opt$htmlPath, sep = "", fill = FALSE, labels = NULL,
                 append = TRUE) {
  if (is.character(file))
    if (file == "")
      file <- stdout()
    else if (substring(file, 1L, 1L) == "|") {
      file <- pipe(substring(file, 2L), "w")
      on.exit(close(file))
    }
    else {
      file <- file(file, ifelse(append, "a", "w"))
      on.exit(close(file))
    }
    .Internal(cat(list(...), file, sep, fill, labels, append))
}

# Function to write code for html head and title
HtmlHead <- function(title) {
  cata("<head>\n")
  cata("<title>", title, "</title>\n")
  cata("</head>\n")
}

# Function to write code for html links
HtmlLink <- function(address, label=address) {
  cata("<a href=\"", address, "\" target=\"_blank\">", label, "</a><br />\n")
}

# Function to write code for html images
HtmlImage <- function(source, label=source, height=500, width=500) {
  cata("<img src=\"", source, "\" alt=\"", label, "\" height=\"", height)
  cata("\" width=\"", width, "\"/>\n")
}

# Function to write code for html list items
ListItem <- function(...) {
  cata("<li>", ..., "</li>\n")
}

TableItem <- function(...) {
  cata("<td>", ..., "</td>\n")
}

TableHeadItem <- function(...) {
  cata("<th>", ..., "</th>\n")
}

# Clear file
#cat("", file=opt$htmlPath)

cata("<html>\n")

cata("<body>\n")
cata("<h3>Heatmap from TCGA data:</h3>\n")


# for (i in 1:nrow(imageData)) {
#   if (grepl("heatmap", imageData$Link[i])) {
#    HtmlImage(imageData$Link[i], imageData$Label[i], width=1000)
#   } else {
#     HtmlImage(imageData$Link[i], imageData$Label[i])
#   }
# }
HtmlImage(imageData$Link, imageData$Label, width=1000)
draw(hm)
#PDFs
# 
# for (i in 1:nrow(linkData)) {
#   if (grepl("heatmap.pdf", linkData$Link[i])) {
#     HtmlLink(linkData$Link[i], linkData$Label[i])
#   }
# }
HtmlLink(linkData$Link, linkData$Label)

cata("<p>Click floppy disc icon associated history item to download ")
cata("all files.</p>\n")



cata("</body>\n")
cata("</html>")




