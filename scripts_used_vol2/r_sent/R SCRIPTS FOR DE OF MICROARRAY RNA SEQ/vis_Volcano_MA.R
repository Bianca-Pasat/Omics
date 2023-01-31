############# VOLCANO PLOTS,MA PLOTS, HEATMAPS, TRAJECTORIES
### FOR RANKPROD


de=RP_TOP
de=rnaseq8
de=rnaseq24
de$ens<-rownames(de)
geneseh<-genes_r8
geneseh<-genes_r24
de<-merge(geneseh, de, by.x=1, by.y="ens")
d0<-de

d0$diffexpressed <- "NO"
d0$diffexpressed[de$FC..class1.class2. > 0.5 & de$P.value < 0.05] <- "UP" ## BECAUSE RP DOES 0 VS 1
d0$diffexpressed[de$FC..class1.class2. < -0.5 & de$P.value < 0.05] <- "DOWN" ## BECAUSE RP DOES 0 VS 1
ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed)) + geom_point() + 
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4")) +
  #geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
  ggtitle("Significantly differentially expressed with FC above 0.5")  + theme_minimal()

de$diffexpressed <- "NO"
de$diffexpressed[de$FC..class1.class2. > 1 & de$P.value < 0.05] <- "UP" ## BECAUSE RP DOES 0 VS 1
de$diffexpressed[de$FC..class1.class2. < -1 & de$P.value < 0.05] <- "DOWN"  ## BECAUSE RP DOES 0 VS 1

library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=de, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=V2)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
  ggtitle("Significantly differentially expressed with FC above 0.5")  + 
  theme_minimal()

ggplot(data=d0, aes(x=FC..class1.class2., y=-log10(P.value), col=diffexpressed, label=hgnc_symbol)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4"))  + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  labs (y = "-log10(P value)", x = "Fold change MKC vs DMSO") +
  ggtitle("Significantly differentially expressed with FC above 0.5")  + 
  theme_minimal()


#### MA PLOTS

de %>% 
  mutate(sig = ifelse(P.value < 0.01, FC..class1.class2., NA)) %>% 
  ggplot(aes( RP.Rsum, FC..class1.class2.)) + ########### CLEARLY THIS IS NOT THE SAME AS BASEMEAN
  geom_point(alpha = 0.1) +
  geom_point(aes(y = sig), colour = "brown", size = 1) +
  scale_x_continuous(trans = "log10") 
#facet_wrap(vars(comparison))


#### CHANGE TO SUIT RANKPROD
# set of candidate genes for clustering
#de$gene_name<-rownames(de)

candidate_genes <- de  %>% 
  filter(P.value < 0.05) %>%    # filter table
  pull(V2) %>%             # extract the gene column as a vector
  unique()       