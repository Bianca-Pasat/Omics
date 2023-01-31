library("sva")

library(ggplot2)
library(tidyverse)

### color palettes
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# read the data
trans_cts <-as.data.frame(y$counts)
trans_cts<-as.data.frame(cleaned_count)
trans_cts<-as.data.frame(adjusted) # ************
trans_cts<-as.data.frame(norm_samples) 
trans_cts<-as.data.frame(adjusted2)
trans_cts<-as.data.frame(adjustedc) 


# isoformSwitchAnalyzer
trans_cts<-as.data.frame(stringtie$counts)
rownames(trans_cts)<-trans_cts[,1]
trans_cts<-trans_cts[,-1]
trans_cts<-trans_cts[,c(4,5,6,10,11,12,1,2,3,7,8,9)]

trans_cts<-as.data.frame(stringtie$abundance)
rownames(trans_cts)<-trans_cts[,1]
trans_cts<-trans_cts[,-1]
trans_cts<-trans_cts[,c(4,5,6,10,11,12,1,2,3,7,8,9)]

sample_info<-sampleInfo
sample_info$replicate<-c(paste("R",1:3,sep=""))

colnames(sample_info)<-c("treatment", "samples","replicate")


## reduced - DM24R3, DM8R1
sample_info<-sampleInfo
sample_info<-sample_info[-c(3,4,9,10),]
sample_info$replicate<-c(rep(c(paste("R",1:2,sep=""),paste("R",2:3,sep="")),2))
colnames(sample_info)<-c("treatment", "sample","replicate")

back<-trans_cts
trans_cts<-back
trans_cts<-trans_cts[,-c(3,6,7,10)]


#exons
#raw
trans_cts<-as.data.frame(norm_samples) 
sample_info<-sampleTable

#raw 24 hours
trans_cts<-as.data.frame(norm_samples24) 
sample_info<-sampleTable24

# clean 8
trans_cts<-as.data.frame(norm_samples8) 
trans_cts<-as.data.frame(cleaned_count) 
sample_info<-sampleTable8

# clean 24
trans_cts<-as.data.frame(cleaned_count) 
sample_info<-sampleTable24

trans_cts<-as.data.frame(adjustedc) 
trans_cts<-as.data.frame(counts) 

#on all
trans_cts<-as.data.frame(cleaned_count) 
sample_info<-sampleTable
trans_cts<-as.data.frame(adjustedc) 



sample_info<-samples_info
#8hours rnaseq batch on all
trans_cts=as.data.frame(adjusted[,c(1,4,5,8,9,10)])
sample_info<-samples_info[c(1,4,5,8,9,10),]

#24 hours rnaseq batch on all
trans_cts=as.data.frame(adjusted[,c(2,3,6,7,11,12)])
sample_info<-samples_info[c(2,3,6,7,11,12),]

#colnames(samples_vector_file)<-c("condition","replicate","sample","treatment")
sample_info <- samples_vector_file

#### 24 hours
# sample_info <- samples_vector_file[c(1:3,7:9),]


# PROTEOMICS
#8 hours
x1=X[c(10:12,22:24),]
trans_cts=as.data.frame(x1)
sample_info<-X_info[c(10:12,22:24),]

############ PCA 

# Create a matrix from our table of counts


pca_matrix <- trans_cts %>% 
  as.matrix() %>% 
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)
# sample_pca <- prcomp(pca_matrix, scale=TRUE)

# Convert matrix to tibble
tibpca<-as_tibble(pca_matrix)



# Convert matrix to tibble - add colnames to a new column called "gene"
tibpca<-as_tibble(pca_matrix, rownames = "sample")

pc_eigenvalues <- sample_pca$sdev^2


pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>% 
mutate(pct = variance/sum(variance)*100) %>% 
mutate(pct_cum = cumsum(pct))


pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct), fill="deeppink4") +
  geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")


# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x


pc_scores <- pc_scores %>% 
as_tibble(rownames = "sample")

# pc_scores %>% 
#   # create the plot
#   ggplot(aes(PC1, PC2, col=sample_info$condition)) +
#   geom_point(size=3) + 
#   geom_text(aes(label = sample), vjust = -1, nudge_y = 1) +
#   labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) + # substitute with % variances
#   ggtitle("PCA. by replicate") +
#   theme_bw() 
  
# plot together

p<-sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC1, y = PC2, col = treatment, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()

p<-sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC1, y = PC2, col = joint, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()

# PCA on all together
b<-sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC1, y = PC2, col = treatment, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","deeppink4","#F0E442", "#D55E00")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()

#### PC2 AND PC3
sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC2, y = PC3, col = treatment, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","deeppink4")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC3", round(pc_eigenvalues$pct[3], digits = 2)), x = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()


#### 3 VARIABLES
a<-sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC1, y = PC2, col = treatment, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()

#4 variables
a<-sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x =PC1, y = PC2, col = treatment, shape = replicate)) +
  scale_colour_manual(values=c("deepskyblue4","cornsilk4","deeppink4","#999999")) +
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) +
  ggtitle("PCA of treatment and replicate") +
  theme_bw()



### MAKE QUICKLY TIBBLES FROM PROCOMB OBJECTS
# library(broom)
# 
# # PC variances (eigen values)
# tidy(sample_pca, matrix = "eigenvalues")
# 
# # variable loadings
# tidy(sample_pca, matrix = "loadings")


# Which genes have the most influence in different pc
pc_loadings <- sample_pca$rotation
pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")


top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  #slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot

library(patchwork)

# Adjust some aspects of each plot
pca_plot <- a + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC scores") +
  theme(legend.position = "none")

loadings_plot <- loadings_plot + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC loadings")

# Put them together
(pca_plot | loadings_plot) + plot_annotation(tag_levels = "A")

# TODO see what happens to pca of data if you remove those particular genes








