###### miRNA analysis
### EdgeR ###
library("edgeR")
library("RColorBrewer")
library("Glimma")
library("limma")
library("biomaRt")
library("sva")
library("pheatmap")

################# Creating a DGEList for use with edgeR

reads_one_per_mirna <- read.csv("~/aitor/microRNA/reads_one_per_mirna.csv") ### load the counts matrix wherever it is
reads_one_per_mirna[,5:28] <- lapply(reads_one_per_mirna[,5:28] , as.numeric)
counts <- reads_one_per_mirna[,5:28]
row.names(counts) <- reads_one_per_mirna$X

samples <- as.data.frame(cbind(rep("a", 24), rep("b", 24), rep("c", 24)))
samples$V1 <- as.factor(c("D8", "D8", "D8", "M8", "M8", "M8", "T8", "T8", "T8", "TM8", "TM8", "TM8", "D24", "D24", "D24", "M24", "M24", "M24", "T24", "T24", "T24", "TM24", "TM24", "TM24"))
samples$V2 <- as.factor(c("R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3", "R1", "R2", "R3"))
samples$V3 <- as.factor(c("8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "8", "24", "24", "24", "24", "24", "24", "24", "24", "24", "24", "24", "24"))
colnames(samples) <- c("condition","replicate", "times")
row.names(samples) <- colnames(counts)

genetable <- reads_one_per_mirna[,2:4]

### create EdgeR object

y <- DGEList(counts=counts, samples=samples, genes=genetable)
names(y)

# COUNT PER MILLION, LOG2 COUNTS
cpm <- cpm(y) 
lcpm <- cpm(y, log=TRUE)

########################## REMOVE LOW EXPRESSED GENES

table(rowSums(y$counts==0)==24)

keep.exprs <- rowSums(y$counts==0) < 12
y <- y[keep.exprs,, keep.lib.sizes=FALSE]

# FILTERED DENSITY PLOT

nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")

plot(density(lcpm[,1]), col=col[1], lwd=2, las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n", y.intersp = 0.5)


lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n", y.intersp = 0.5)


######################## Normalisation by the method of trimmed mean of M-values (TMM)

wee <- log2(y$counts)
boxplot(wee, las=2, col=col, main="")
title(main="Log2 Raw data",ylab="Log-cpm")

lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Log2 CPM data",ylab="Log-cpm")

y <- calcNormFactors(y)
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

# MDS plots

plotMDS(lcpm, top = 1000, labels = NULL, col = as.numeric(y$samples$condition), cex = 2)
title(main="Treatment groups")
plotMDS(lcpm, top = 1000, labels = NULL, col = as.numeric(y$samples$replicate), cex = 2)
title(main="Replicate groups")

glMDSPlot(lcpm, labels=paste(condition, replicate, sep="_"), groups=y$samples[,c(2,5)], launch=FALSE) # Glimma

# sample-to-sample distances

sampleDists <- dist(t(lcpm))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(lcpm)
colnames(sampleDistMatrix) <- colnames(lcpm)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

############################### BATCH CORRECTION. SVA estimation of surrogate variables

mod1 <- model.matrix(~0 + condition, y$samples)
mod0 <- model.matrix(~1, y$samples)

svobj <- svaseq(cpm(y), mod1, mod0) 
design <- cbind(mod1, svobj$sv)

### "Clean" gene expression data, with the calculated sva 
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- cleanY(cpm(y), mod1, svobj$sv)
log_cleaned_count <- log2(cleaned_count)


# Re-Do PCA, hierarchical clust, etc.
plotMDS(log_cleaned_count, top = 1000, labels = NULL, col = as.numeric(y$samples$condition), cex = 2)
title(main="Treatment groups")
plotMDS(log_cleaned_count, top = 1000, labels = NULL, col = as.numeric(y$samples$replicate), cex = 2)
title(main="Replicate groups")

pheatmap(log_cleaned_count, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, scale = "row" , annotation_col=samples)

############################ DIFFERENTIAL EXPRESSION ANALYSIS EdgeR
#each gene will get its own unique dispersion estimate, but the common dispersion is still used in the calculation.

y <- estimateDisp(y, design)

# mean-variance plot. raw variances of the counts (grey dots), the variances using the tagwise

meanVarPlot <- plotMeanVar( y , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

### DEG
fit <- glmFit(y, design)
D8vsM8 <- glmLRT(fit, contrast = c(0,-1,0,1,0,0,0,0,0,0,0,0)) # -1 is the denominator,, as many 
D8vsM8_DEG <- topTags(D8vsM8, n=nrow(y), p.value=0.1)
res_D8vsM8 <- topTags(D8vsM8, n=nrow(y))
resultados_D8vsM8 <- as.data.frame(res_D8vsM8)
topTags(D8vsM8) # just the top 10 by default

table(resultados_D8vsM8$FDR < 0.1)

results_edgeR_8 <- Reduce(function(x, y) merge(x, y, by= "gene.id"), list(resultados_D8vsM8, resultados_T8vsTM8, resultados_D8vsT8))

### Visualizing results
# Export to Glimma

glMDPlot(D8vsM8, 
         counts=y$counts, 
         anno=y$genes, 
         groups=y$samples$condition,
         status= probando,
         samples=colnames(y))

# Abundance of top DEG genes

hist( resultados_D8vsM8[1:100,"logCPM"] , breaks=10 , xlab="Log Concentration" ,
      col="red" , freq=FALSE , main="Poisson: Top 100" )

# MA plot

plotSmear(D8vsM8, de.tags=D8vsM8_DEG$table$gene.id, cex = 1, col = "grey10" )

### Volcano plot 

resultados_D8vsM8$TREND = "N.C" 
resultados_D8vsM8$TREND[resultados_D8vsM8$logFC > 0.25 & resultados_D8vsM8$FDR < 0.1 ] <- "UP"   ### adjust resultados_D8vsM8$logFC with whatever log2FC cutoff you want, ex=0.25
resultados_D8vsM8$TREND[resultados_D8vsM8$logFC < -0.25 & resultados_D8vsM8$FDR < 0.1 ] <- "DOWN" 

TFs_8h_volc <- ggplot(resultados_D8vsM8, aes(resultados_D8vsM8$logFC, -log10(resultados_D8vsM8$FDR)))
TFs_8h_volc + geom_point(aes(colour = resultados_D8vsM8$TREND, size = 1, alpha = 0.5)) + labs(x = "log2 Fold Change", y = "-log10 p.value", element_text(face = "bold", angle = 0)) + scale_colour_manual(values = c("green", "black", "red")) + geom_vline(xintercept = 0.25 , color = "brown", linetype = 2) +  geom_vline(xintercept = -0.25 , color = "brown", linetype = 2) + geom_hline(yintercept = 1 ,linetype = 2) + geom_hline(yintercept = 0 ) + geom_vline(xintercept = 0)+ theme(panel.background = element_rect(fill = "white", colour = "grey50"), panel.grid.major = element_line(colour = "gray85"))


