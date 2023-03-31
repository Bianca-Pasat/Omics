#Rscript load_heatmap.R -c '/Users/user/GDCdata/TCGA-COAD.RData' -g '/Users/user/NUIG-HCT116/wptb_comon_network/Gene_Prioritization.tsv' -o F -html /Users/user/testBIANCA.html
rm(list = ls()) 

suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(getopt))

args = commandArgs(trailingOnly=TRUE)


#cols=eval(parse(text=args[1]))
# cancer<-as.character(args[1])
# genesign<-read.table(args[2], header = T,  sep="\t")
# htmlPath<-args[3]
# outPath<-args[4]


# Initialise data for html links and images, data frame with columns Label and
# Link
linkData <- data.frame(Label=character(), Link=character(),
                       stringsAsFactors=FALSE)
imageData <- data.frame(Label=character(), Link=character(),
                        stringsAsFactors=FALSE)



# Get options, using the spec as defined by the enclosed list.
# Read the options from the default: commandArgs(TRUE).

spec <- matrix(c(
  "cancer", "c", 1, "character",
  "genesign", "g", 1, "character",
  "optional_parameters", "o", 1,"logical",
  "cluster_rows", "cr",2,"logical",
  "cluster_columns", "cc",2,"logical",
  "show_row_names", "srn",2,"logical",
  "show_column_names", "scn",2,"logical",
  "show_row_dend", "srd",2,"logical",
  "show_column_dend", "scd",2,"logical",
  "row_names_max_width", "rnmw", 2, "double",
  "clustering_distance_rows", "cdr", 2, "character",
  "clustering_method_rows", "cmr", 2, "character",
  "clustering_distance_columns", "cdc", 2, "character",
  "clustering_method_columns", "cmc", 2, "character",
  "row_dend_width", "rdw", 2, "double",
  "row_km", "rkm", 2, "integer",
  "column_km", "ckm", 2, "integer",
  "km", "km", 2, "integer",
  "heatmap_width", "w", 2, "integer",
  "heatmap_height", "h", 2, "integer",
  "survival", "s", 1, "character",
  "html","html","2","character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)
#  

#genesign<-read.table("/Users/user/NUIG-HCT116/wptb_comon_network/Gene_Prioritization.tsv",  header = T,  sep="\t")
load(opt$cancer)
genesign<-read.table(opt$genesign, header=T,row.names=1,sep="\t")

genesign <- as.character(rownames(genesign))


genesign.FC <- as.matrix(genesign[,1])
rownames(genesign.FC) <- rownames(genesign)

#keep the genes of signature that are present in exp.filt for further analysis
genesign <-  genesign[genesign %in% rownames(exp.filt)]
genesign.FC <- genesign.FC[rownames(genesign.FC) %in% genesign, ]

#map the genesign on logCPM.matrix
logCPM.genesign <- logCPM[genesign.FC, ]

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

TCGAanalyze_survival(data = colData(exp),
                     clusterCol = clusterCol,
                     #filename = "zsurvival.png",#or "survivalPlot.png"
                     #height = 10,
                     #width = 10,
                     dpi = 1200,
                     legend = "Survival Plot",
                     main = "Kaplan-Meier Overall Survival Curves",
                     pvalue = TRUE,
                     risk.table = TRUE,
                     conf.int= FALSE)
linkData <- rbind(linkData, "survival.pdf")
imageData<- rbind(imageData, "survival.pdf")
dev.off()



###heatmap
colData(exp)$heat <- as.numeric("0")
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"heat"] <- as.numeric("-1")
colData(exp)[score >= high.cut,"heat"] <- as.numeric("1")


testx<-matrix(colData(exp)$heat, nrow=nrow(exp), ncol=length(colData(exp)$barcode))
rownames(testx) <- rownames(exp)
colnames(testx) <- colnames(logCPM.genesignx)
testx<-testx[genesignx,]

test<-testx
##annotation and parameters for heatmap

col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5
clinical_patient_Cancer<-colData(se)
clinical_patient_Cancer<-clinical_patient_Cancer[rownames(clinical_patient_Cancer) %in% colnames(test),]


signa<-data.frame(rownames(test))

col<-HeatmapAnnotation(tissue_type=clinical_patient_Cancer[,4])

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")


if ((opt$optional_parameters)){
  cluster_rows=opt$cluster_rows
  cluster_column=opt$cluster_column
  show_row_name=opt$show_row_name
  show_column_name=opt$show_column_name
  show_row_dend=opt$show_row_dend
  show_column_dend=opt$show_column_dend
  row_names_max_width=opt$row_names_max_width
  clustering_distance_rows=opt$clustering_distance_rows
  clustering_method_rows=opt$clustering_method_rows
  clustering_distance_columns=opt$clustering_distance_columns
  clustering_method_column=opt$clustering_method_column
  row_dend_width=opt$row_dend_width
  row_km = opt$row_km
  column_km = opt$column_km
  km=opt$km
  heatmap_width=opt$heatmap_width
  heatmap_height=opt$heatmap_height
}else{
  cluster_rows=T
  cluster_column=T
  show_row_name=T
  show_column_name=F
  show_row_dend=FALSE
  show_column_dend=TRUE
  row_names_max_width=0.1
  clustering_distance_rows="euclidean"
  clustering_method_rows="ward.D"
  clustering_distance_columns="euclidean"
  clustering_method_column="ward.D"
  row_dend_width=10
  row_km=1
  column_km=1
  km=1
  heatmap_width=15
  heatmap_height=15
}


#heatmap
pdf("heatOutPdf.pdf")
Heatmap(test, name="expression",
        col=col_fun,
        cluster_rows=cluster_rows,
        cluster_columns = cluster_column,
        row_names_side = "left",
        show_row_names = show_row_name,
        #column_order = order(clinical_patient_Cancer$shortLetterCode),
        #row_order = order(rownames(logCPM.genesign)),
        show_column_names = show_column_name,
        #column_names_gp=gpar(cex=fontsize),
        show_row_dend = show_row_dend,
        show_column_dend=show_column_dend,
        row_names_gp=gpar(cex=fontsize+row_names_max_width),
        row_names_max_width = unit(row_names_max_width, "cm"),
        clustering_distance_rows =clustering_distance_rows,
        clustering_method_rows = clustering_method_rows,
        clustering_distance_columns =  clustering_distance_columns,
        clustering_method_columns = clustering_method_column,
        row_dend_width = unit(row_dend_width, "mm"),
        row_km = row_km,
        column_km = column_km,
        km=km,
        left_annotation=row,
        bottom_annotation = col,
        heatmap_width = unit(heatmap_width, "cm"),
        width = NULL,
        heatmap_height = unit(heatmap_height, "cm"))
# draw(hm)
# linkName <- paste0("heatmap", ".pdf")
# linkAddr <- paste0("heatmap", ".pdf")
linkData <- rbind(linkData, "heatOutPdf.pdf")
imageData<- rbind(imageData, "heatOutPdf.pdf")
invisible(dev.off())




################################################################################
### HTML Generation
################################################################################
# Create cata function: default path set, default seperator empty and appending
# true by default (Ripped straight from the cat function with altered argument
# defaults)


cata <- function(..., file = opt$html, sep = "", fill = FALSE, labels = NULL,
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

# # Function to write code for html head and title
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
cat("", file=opt$html)

cata("<html>\n")

cata("<body>\n")

cata("<h3>Plots from TCGA data:</h3>\n")

# 
for (i in 1:nrow(imageData)) {
  HtmlImage(imageData$Link[i], imageData$Label[i], width=1000)
}

# #PDFs
# # 
for (i in 1:nrow(linkData)) {
  HtmlLink(linkData$Link[i], linkData$Label[i])
}

# 
cata("<p>Click floppy disc icon associated history item to download ")
cata("all files.</p>\n")



cata("</body>\n")
cata("</html>")
