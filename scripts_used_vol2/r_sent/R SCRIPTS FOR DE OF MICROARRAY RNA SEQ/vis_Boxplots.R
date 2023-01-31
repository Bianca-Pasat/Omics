trans_cts <- as.data.frame(y$counts)
trans_cts<- as.data.frame(norm_samples)
sample_info<-y$samples

##isoforms
trans_cts<- as.data.frame(stringtie$abundance)
rownames(trans_cts)<-trans_cts[,1]
trans_cts<-trans_cts[,-1]

trans_cts<- as.data.frame(stringtie$counts)
rownames(trans_cts)<-trans_cts[,1]
trans_cts<-trans_cts[,-1]

sample_info<-sampleInfo
colnames(sample_info)<-c("condition","sample")
sample_info$times<-rep(c(rep(8,3),rep(24,3)),2)
sample_info$replicate<-rep(paste("R",1:3,sep=""),4)

###exons
trans_cts<- as.data.frame(counts1)
trans_cts<- as.data.frame(counts)

trans_cts<- as.data.frame(norm_samples)
sample_info<-sampleTable

trans_cts<- as.data.frame(norm_samples8)
sample_info<-sampleTable8

#after correction
trans_cts<- as.data.frame(cleaned_count)
sample_info<-sampleTable24
sample_info<-sampleTable8

# PROTEOMICS
trans_cts<- as.data.frame(norm_proteins)
trans_cts<-as.data.frame(proteins)
sample_info<-sample_info_proteins
colnames(sample_info)=c("whole_name","hours", "condition","replicate","joint","sample")

###microarray
trans_cts <- as.data.frame(y$counts)
trans_cts<- as.data.frame(exp)
sample_info<-y$samples

##parameterization
vis1<-1:6

########## Boxplots 
raw_cts_long <- trans_cts %>% 
  pivot_longer(DMSO1:MKC6, names_to = "sample", values_to = "expression")

raw_cts_long <- trans_cts %>% 
  pivot_longer(1:12, names_to = "sample", values_to = "expression")

raw_cts_long <- trans_cts %>% 
  pivot_longer(1:24, names_to = "sample", values_to = "expression")

raw_cts_long <- trans_cts %>% 
  pivot_longer(1:6, names_to = "sample", values_to = "expression")


# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))


########### ADD AFTER NORMALIZATION


raw_cts_long <- trans_cts %>% 
  pivot_longer(DMSO_8_1:MKC8866_24_1, names_to = "sample", values_to = "expression")

raw_cts_long <- trans_cts %>% 
  pivot_longer(D8R1:M24R3, names_to = "sample", values_to = "expression")

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))


## whole once
raw_cts_long <- trans_cts %>% 
  pivot_longer(1:12, names_to = "sample", values_to = "expression")

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Make a boxplot
raw_cts_long %>%
  ggplot(aes(as.character(condition), log2(expression),fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of raw counts before normalization") +
  theme_bw()

raw_cts_long %>%
  ggplot(aes(as.character(hours), log2(expression),fill = condition)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of counts after normalization") +
  theme_bw()

boxplotidi<-function(){
  raw_cts_long <- trans_cts %>% 
    pivot_longer(vect, names_to = sample, values_to = "expression")
  raw_cts_long <- full_join(raw_cts_long, sample_info, by = (sample))
  raw_cts_long %>%
    ggplot(aes(as.character(times), log2(expression),fill = treatment)) + 
    geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
    facet_grid(cols = vars(replicate)) +
    labs (y = "log2 expression", x = "time points") +
    ggtitle("Boxplot of raw counts before normalization") +
    theme_bw()
}

raw_cts_long %>%
  ggplot(aes(as.character(replicate), log2(expression),fill = condition)) + 
  geom_boxplot() + scale_fill_manual(values=cbp2) + 
  facet_grid(cols = vars(times)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of raw counts before normalization") +
  theme_bw()


#exons after correction
raw_cts_long %>%
  ggplot(aes(as.character(treatment), log2(expression),fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of raw counts after batch correction") +
  theme_bw()


#exons
raw_cts_long %>%
  ggplot(aes(as.character(times), log2(expression),fill = condition)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of counts after normalization") +
  theme_bw()

raw_cts_long %>%
  ggplot(aes(as.character(times), log2(expression),fill = condition)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of log2 counts before normalization") +
  theme_bw()


raw_cts_long %>%
  ggplot(aes(as.character(times), log2(expression),fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "time points") +
  ggtitle("Boxplot of log2 counts after normalization") +
  theme_bw()


raw_cts_long %>%
  ggplot(aes(expression, log2(expression),fill = treatment)) + 
  geom_boxplot()  + 
  facet_grid(cols = vars(replicate)) +
  ggtitle("Boxplot of log2 counts after normalization") +
  theme_bw()

raw_cts_long %>%
  ggplot(aes(expression,expression,fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "expression") +
  ggtitle("Boxplot of log2 counts after normalization with Upper Quartile") +
  theme_bw()

raw_cts_long %>%
  ggplot(aes(expression,expression,fill = treatment)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(replicate)) +
  labs (y = "log2 expression", x = "expression") +
  ggtitle("Boxplot of log2 counts after normalization with Upper Quartile") +
  theme_bw()



########### ADD cleaning
raw_cts_long <- trans_cts %>% 
  pivot_longer(DMSO1:MKC3, names_to = "sample", values_to = "expression")

# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))

# Make a boxplot
raw_cts_long %>%
  # make sure minute is specified as a factor
  ggplot(aes(expression, log2(expression),fill = treatment)) + 
  geom_boxplot() + scale_fill_manual(values=c("deepskyblue4","deeppink4")) + 
  facet_grid(cols = vars(replicate)) +
  ggtitle("Boxplot of log2 counts after cleaning") +
  theme_bw()

