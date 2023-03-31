# load packages
library(tidyverse)


### color palletes
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# read the data
trans_cts <-as.data.frame(y$counts)

colnames(samples_vector_file)<-c("condition","replicate","sample","treatment")
sample_info <- samples_vector_file


raw_cts_long <- trans_cts %>% 
  pivot_longer(DMSO1:MKC3, names_to = "sample", values_to = "expression")

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Make a boxplot
raw_cts_long %>%
  # make sure minute is specified as a factor
  ggplot(aes(expression, log10(expression),fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=cbp2) + 
  facet_grid(cols = vars(replicate))

# Create a matrix from our table of counts
pca_matrix <- trans_cts %>% 
  # make the "gene" column become the rownames of the table
  #column_to_rownames("gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)

# Convert matrix to tibble
tibpca<-as_tibble(pca_matrix)


# Convert matrix to tibble - add colnames to a new column called "gene"
tibpca<-as_tibble(pca_matrix, rownames = "sample")

pc_eigenvalues <- sample_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>% 
# add a new column with the percent variance
mutate(pct = variance/sum(variance)*100) %>% 
# add another column with the cumulative variance explained
mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct), fill="deeppink4") +
  geom_line(aes(y = pct_cum, group = 1), col="deepskyblue4") + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")


# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x


pc_scores <- pc_scores %>% 
# convert to a tibble retaining the sample names as a new column
as_tibble(rownames = "sample")

# print the result
pc_scores

pc_scores %>% 
  # create the plot
  ggplot(aes(PC1, PC2, col=sample_info$pca1)) +
  geom_point(size=3) + 
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) +
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) + # substitute with % variances
  ggtitle("PCA. by replicate") +
  theme_bw() 

# plot together
sample_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  full_join(sample_info, by = "sample") %>% 
  # make the plot
  ggplot(aes(x = PC1, y = PC2, 
             col = treatment, shape = replicate)) + 
  geom_point(size=5) +
  geom_text(aes(label = sample), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(pc_eigenvalues$pct[2], digits = 2)), x = paste("PC1", round(pc_eigenvalues$pct[1], digits = 2))) 


# Which genes have the most influence in different pc
pc_loadings <- sample_pca$rotation
pc_loadings <- pc_loadings %>% 
as_tibble(rownames = "gene")

# print the result
pc_loadings

top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC3) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes

# TODO see what happens to pca of data if you remove those particular genes


### MAKE QUICKLY TIBBLES FROM PROCOMB OBJECTS
# library(broom)
# 
# # PC variances (eigen values)
# tidy(sample_pca, matrix = "eigenvalues")
# 
# # variable loadings
# tidy(sample_pca, matrix = "loadings")

