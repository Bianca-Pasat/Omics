######### Clustering options

# Summarise counts
trans_cts<-as.data.frame(adjusted)
trans_cts$gene<-rownames()
trans_cts_mean <- trans_cts %>% 
  # convert to long format
  pivot_longer(cols = DMSO1:MKC_24_3, names_to = "sample", values_to = "expression")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>% 
  # filter to retain only genes of interest
  filter(gene %in% candidate_genes) %>% 
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(expression_scaled = (expression- mean(expression))/sd(expression)) %>% 
  # for each gene, sample and treatment
  group_by(gene,sample,condition) %>%
  # calculate the mean (scaled) cts
  summarise(mean_expression_scaled = mean(expression_scaled), 
            nrep = 3) %>% 
  ungroup()


# Create a matrix
hclust_matrix <- trans_cts %>% 
  select(-gene) %>% 
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene

hclust_matrix <- hclust_matrix[rownames(hclust_matrix) %in% candidate_genes, ]

hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


gene_dist <- dist(hclust_matrix)

gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 10, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 5) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
