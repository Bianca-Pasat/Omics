library("depmap")
library("ExperimentHub")
library(tidyverse)

## create ExperimentHub query object
eh <- ExperimentHub()
## snapshotDate(): 2022-10-24
qu<-query(eh, "depmap")

metadata <- eh[["EH2266"]]
pancs<-c("ACH-000107","ACH-000060","ACH-000138","ACH-000094","ACH-000155","ACH-000535","ACH-000222","ACH-000205","ACH-000601","ACH-000265","ACH-000085","ACH-000307")
si_pancs<-

clines<-depmap_TPM()
panclines_meta<-clines[clines$depmap_id %in% pancs,]
panclines1<-clines[clines$depmap_id %in% pancs,][,c(1,3,5)]
panclines<-as.data.frame(panclines1  %>% pivot_wider(
                            names_from = depmap_id ,
                            values_from = rna_expression))
rm(panclines1)

rownames(panclines)<-panclines$gene_name
panclines<-panclines[,-1]
write.table(panclines_meta,"/home/ben/Desktop/pancreatic_cell_lines_metadata.txt",sep="\t",quote=F, row.names=T)
write.table(panclines,"/home/ben/Desktop/pancreatic_cell_lines.txt",sep="\t",quote=F, row.names=T)

si_pancs<-as.data.frame(unique.data.frame(panclines_meta[,c(1,6)]))
si_pancs$condition<-as.numeric(rep(1, length(si_pancs$depmap_id)))
rownames(si_pancs)<-si_pancs$depmap_id
# design<-model.matrix(~0+condition, si_pancs)
# 
# ### filter and logcpm/normalize
# preproc<-function(counts,si,design,rpnorm){
#   ### make object, filter and normalize
#   counts <- mutate_all(counts, function(x) as.numeric(as.character(x)))
#   y <- DGEList(counts=counts, samples=si)
#   # REMOVE LOW EXPRESSED GENES
#   keep.exprs <- filterByExpr(y,design=design)
#   y <- y[keep.exprs,, keep.lib.sizes=FALSE]
#   ####### normalization
#   if(rpnorm){
#     y$counts<-function.norm(y$counts)
#   }
#   else if(!rpnorm){
#     y <- calcNormFactors(y, method = "upperquartile")
#   }
#   return(y)
# }
# logcpmPanc<-preproc(panclines,si_pancs,design, FALSE)[[1]]

#create DGEList object (edgeR)
y.exp <- DGEList(counts = panclines)

#TMM normalization
y.exp <- calcNormFactors(y.exp, method = "TMM")

###LOGCPM
logCPM <- cpm(y.exp, prior.count=5, log=TRUE)
logCPM  <- TCGAanalyze_Filtering(tabDF = logCPM,
                                 method = "quantile",qnt.cut =  0.25)



### load gene signature
genesign<-read.table("/home/ben/Desktop/xbp1_survival/cytokines_paad.txt", sep="\t")
genesign<-genesign$V1
# genesign<-read.table("/home/ben/Desktop/xbp1_survival/paad_cytokines_GenePriori.tsv", sep="\t")
# genesign<-genesign$V2
### put FC =1 because you don't know if cytokines going down are ridd targets.
genesign.FC<-c(rep(1,length(genesign)))
#### calculate with KV method
logCPM.genesign <- logCPM[rownames(logCPM) %in% genesign, ]

# logCPM.genesignr <- logCPM[genesignr, ]


#calculate the 25th, 50th and 75th quantile per row
probs <- c(0.25,0.5,0.75)
# probs <- c(0.33,0.5,0.67)
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

score <- colSums(aux)/length(genesign)

#set a threshold to group patients  e.g.0.25 
threshold <- 0.25
threshold <- 0.3

#define the groups for comparison based on the user-defined threshold
colData(exp)$level.group.xbp1 <- "Mid expression"
min.cut <- max(sort(score)[1:(length(score) * threshold)])
high.cut <- min(sort(score, decreasing = T)[1:(length(score) * threshold)])
colData(exp)[score <= min.cut,"level.group.xbp1"] <- "Low expression"
colData(exp)[score >= high.cut,"level.group.xbp1"] <- "High expression"








test<-aux


col_fun = colorRamp2(seq(min(test), max(test), length = 3), c("blue", "white", "red"))
fontsize=0.5
signa<-data.frame(rownames(test))
col<-HeatmapAnnotation(cell_lines=si_pancs$cell_line)

row<-HeatmapAnnotation(signature=signa[,1] ,which="row")
hm<-Heatmap(test, name="expression",
            col=col_fun,
            cluster_rows=T,
            cluster_columns = T,
            #row_names_side = "left",
            show_row_names = T,
            #column_order = order(clinical_patient_Cancer$shortLetterCode),
            #row_order = order(rownames(logCPM.genesign)),
            show_column_names = F,
            #column_names_gp=gpar(cex=fontsize),
            show_row_dend = TRUE,
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
            #left_annotation=row,
            bottom_annotation = col,
            heatmap_width = unit(20, "cm"),
            width = NULL,
            heatmap_height = unit(20, "cm"))
draw(hm)
