############## trends of genes
trans_cts_cluster <- trans_cts_mean %>% 
  inner_join(gene_cluster, by = "gene")

trans_cts_cluster %>% 
  arrange(gene) %>%
  ggplot(aes(condition,mean_expression_scaled,  group=1)) +
  geom_line(aes(group= gene)) +
  #facet_grid(rows = vars(sample), cols = vars(cluster)) +
  theme_bw()

trans_cts_cluster %>% 
  ggplot(aes(mean_expression_scaled)) +
  geom_line(aes(group = gene)) +
  facet_grid(rows = vars(treatment), cols = vars(cluster)) +
  theme_bw()

trans_cts_cluster %>% 
  ggplot(aes(minute, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  facet_grid(rows = vars(strain), cols = vars(cluster))

