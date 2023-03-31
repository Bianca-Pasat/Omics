#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
################# CORRELATION MATRICES HAVE ROWNAMES GENES AND MIRNAS AND COLNAMES GENES AND MIRNAS AND VALUES ARE THE CORRELATION 

cor_0<-readRDS("/home/ben/Desktop/mirna_bianca/correlation_8_rna")

cor_0<-as.data.frame(corr_0[,c(1,2,3)])
cor_01<-as.data.frame(corr_0[,c(2,1,3)])
colnames(cor_0)=c("mirna","mrna","correlation")
colnames(cor_01)=c("mirna","mrna","correlation")
cor_0<-rbind(cor_0, cor_01)

cor_0$correlation<-as.numeric(cor_0$correlation)
hey<-cor_0 %>% pivot_wider(names_from = mrna, values_from = correlation)
hey<-cococo %>% pivot_wider(names_from = mrna, values_from = correlation)
hey[is.na(hey)]<-0
input.matrix = as.matrix(read.table(args[1]))
input.matrix = as.matrix(hey)
class(input.matrix)="numeric"

rownames(input.matrix)=hey$mirna
#colnames(input.matrix)=input.matrix[1,]
input.matrix=input.matrix[,-1]
input.matrix=input.matrix[-1,-1]

dm8<-dem8[colnames(hey),]
dmi8<-dedemi8[hey$mirna,]
dm8<-as.data.frame(dm8[,-c(1,3)], row.names = rownames(dm8))
dmi8<-as.data.frame(dmi8[,-c(1,3)], row.names = rownames(dmi8))
colnames(dm8)="fc"
colnames(dmi8)="fc"
cormat<-as.data.frame(rbind(dm8, dmi8))

input.matrix=as.matrix(cormat)
################# CORRELATION MATRICES HAVE ROWNAMES GENES AND MIRNAS AND COLNAMES GENES AND MIRNAS AND VALUES ARE THE CORRELATION 

########### correlation matrix
# cor_pearson<-cor(fused_sigpatients, method = c("pearson"))
# cor_kendall<-cor(fused_sigpatients, method = c("kendall"))
# cor_spearman<-cor(fused_sigpatients, method = c("spearman"))

library('Hmisc')
# TYPE OF CORRELATION: "pearson"(default), "spearman"  
# Missing values are deleted in pairs (case-wise deletion)
typeof.correlation = as.character("pearson")
results <- rcorr(input.matrix, type = typeof.correlation)
results <- rcorr(as.matrix(input.matrix), type = typeof.correlation)
results <- rcorr(input.matrix[,c(1:13)], type = typeof.correlation)
# Extract the correlation coefficients
write.table(results$r, "matrix_of_correlations.tsv", sep = "\t")
# Extract p-values
write.table(results$P, "matrix_of_Pvalues.tsv", sep = "\t")

library('corrplot')
# ORDERING OF CORR MATRIX BASED ON : "FPC", "hclust" or "original" (default)
# DISPLAY TYPE: "full" (default), "upper" (upper triangular matrix)
# SIGNIFICANT LEVEL: 0.01 (default, numeric value) 
# if the p-value in p-mat is bigger than sig.level, then the corresponding correlation 
# coefficient is regarded as insignificant and left blank
corr.display = as.character(args[3])
corr.display = as.character("upper")
corr.order = as.character(args[4])
corr.order = as.character("FPC")
corr.siglevel = as.numeric(args[5])
corr.siglevel = as.numeric(0.05)

corr.signif = as.logical(args[6])
corr.signif = as.logical(TRUE)
if (!corr.signif) {
  corr.insig = "n"
} else {
  corr.insig = "blank"
}

png(height=750, width=700, file="correlation_matrix.png", type = "cairo")
corrplot(results$r, type="lower", order=corr.order, p.mat = results$P, 
         sig.level = corr.siglevel, insig = corr.insig,
         tl.cex = 0.7, tl.srt = 30, tl.col = "black", method = "square")
corrplot(results$r, type="lower", order="original", p.mat = results$P, 
         sig.level = corr.siglevel, insig = corr.insig,
         tl.cex = 0.7, tl.srt = 30, tl.col = "black", method = "square")
dev.off()

library(pheatmap)
col <- colorRampPalette(c("red", "white", "blue"))(20)
breaks <- seq(1, -1, length.out = 21)

pheatmap(results$r, show_rownames = F, show_colnames = F, col = col,
         cellheight = 3, cellwidth = 3, breaks = breaks,
         filename = "corrcoef_heatmap.png", type = "cairo")


sig.results <- matrix( 0, dim(results$r)[1], dim(results$r)[2])
for (i in 1:nrow(results$r)) {
  for(j in 1:ncol(results$r)) {
    sig.results[i,j][results$P[i,j]<=corr.siglevel] <- results$r[i,j]
}}
write.table(sig.results, "matrix_of_significant_correlations.tsv", sep = "\t")
pheatmap(sig.results, show_rownames = F, show_colnames = F, col = col,
         cellheight = 3, cellwidth = 3, breaks = breaks,
         filename = "significant_corrcoeff_heatmap.png", type = "cairo")

######################## PCA
# Compute PCA by performing SVD of the centered and scaled data matrix
pca <- prcomp(input.matrix, scale. = TRUE) #center = TRUE (default)
write.table(data.frame(summary(pca)$importance), "PCA_summary.tsv", sep="\t")
write.table(data.frame(pca$rotation), "eigenvectors_matrix.tsv", sep="\t")

library(factoextra)
# Visualize the eigenvalues of dimensions (amount of the variation explained by each PC)
fviz_eig(pca, addlabels = TRUE) # main="..."
ggsave(filename = "scree_plot.png", type = "cairo")

# Graph of individuals. Individuals with a similar profile are grouped together.
fviz_pca_ind(pca, repel = TRUE, # Avoid text overlapping
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("blue","yellow","red")) # else automatic coloring by cos2
ggsave(filename = "individuals_PCA.png", type = "cairo")

# Graph of variables. 
# Positive correlated variables point to the same side of the plot. 
# Negative correlated variables point to opposite sides of the graph.

var.selection = as.logical(args[8])
if (!var.selection) { 
  var.contrib = NULL
} else { 
  var.contrib = as.integer(args[9])
}

fviz_pca_var(pca, repel = TRUE,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("blue", "yellow", "red"),
             select.var = list(contrib = var.contrib))
ggsave(filename = "variables_PCA.png", type = "cairo")

# Biplot of individuals and variables
biplot.label = as.character(args[7]) # all, ind, var
fviz_pca_biplot(pca, repel = TRUE, label = biplot.label)
ggsave(filename = "PCA_biplot.png", type = "cairo")

# Eigenvalues
eig.val <- get_eigenvalue(pca)
write.table(data.frame(eig.val), "eigenvalues_matrix.tsv", sep="\t")

