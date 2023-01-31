x <- scan("/home/bianca/Downloads/ml_covid/nextclade_after_june_small_tab.tsv", what = "", sep = "\n")
x <- strsplit(x, "[ \t]+") # split string by white space
max.col <- max(sapply(x, length))

## option 1
## specify col.names as ?read.table suggests
cn <- paste("V", 1:max.col, sep = "")
after_june <- read.table("/home/bianca/Downloads/ml_covid/nextclade_after_june_small_tab.tsv", fill = TRUE, col.names = cn)
after_june <- read.table("/home/bianca/Downloads/ml_covid/nextclade_after_june_small_tab.tsv", fill = TRUE, col.names = sites)
#rownames(after_june)=after_june$V1
rownames(after_june)=after_june$V1
after_june1<-replace(after_june,"","NA")
long_after_june<-after_june %>%
  pivot_longer(!V1,names_to = "idk", values_to = "mutations",values_drop_na = TRUE)

df1<-c()
df1<-data.frame(nrows=1)
for (i in 1:829){
  df1<-append(df1, after_june[i,c(2:56)])
}
sites<-unlist(dfone)
sites=sites[-c(1,2)]

ex<-read.table("/home/bianca/Downloads/ml_covid/inputs_for_nextclade/after_june_numcer/covid_annot.csv",sep=",")
colnames(ex)=ex[1,]
ex=ex[-1,]
library(tidyverse)
#Dataframe

df <- data.frame(variant=unique(ex$variant),stringsAsFactors = F)

#Code
ex1 <- ex %>% select(1,8) %>% mutate(Value=1) %>%
  full_join(df) %>%
  fill(sample) %>%
  pivot_wider(names_from = variant,values_from=Value) %>%
  replace(NULL,0)

ex <- data.frame(billno = c(715851, 715851, 715851,715852, 715852, 715852, 715852, 715852, 715852), signatories = c("Ben", "Lisa", "Roger", "Louise", "Macy", "John", "Jake", "James", "Ben"))

Senatornames <- c("Ben", "Lisa", "Roger", "Louise", "Macy", "John", "Jake", "James", "Julian", "Ayn")


library(tidyverse)
#Dataframe
df <- data.frame(signatories=Senatornames,stringsAsFactors = F)
#Code
ex1 <- ex %>% mutate(Value=1) %>%
  full_join(df) %>%
  fill(billno) %>%
  pivot_wider(names_from = signatories,values_from=Value) %>%
  replace(is.na(.),0)