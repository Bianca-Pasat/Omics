
#####3IF YOU WANT TO MAKE POSITIVE VALUES
adjustedcPos<-adjustedc - (min(adjustedc))
svaob<-adjustedcPos


svaob<-y$counts

svaob<-y$counts[,(1:6)]
svaob<-adjustedc
svaob<-norm_samples

#exon usage

svaob<-as.matrix(counts)
countb<-counts
#combat and sva
svaob<-adjustedc - min(adjustedc)



#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod0 <- model.matrix(~1, y$samples)

mod0 <- model.matrix(~1, samples_vector_file)


#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + treatment, y$samples)


#exon usage
mod0 <- model.matrix(~1, sampleTable)

#exon
mod1 <- model.matrix(~0 + condition, sampleTable)

############ IF YOU ARE TAKING THEM ALL TOGETHER THEN YOU NEED TO CHANGE THE TREATMENT
mod1 <- model.matrix(~0 + treatment+times, y$samples)
mod1 <- model.matrix(~0 + treat, y$samples)
mod1 <- model.matrix(~0 + treatment, samples_vector_file)

# ###### alternative
# 
# ### NULL MODEL ONLY AGJUSTMENT VARIABLES (not the ones that we are interested in) 
# mod0 <- model.matrix(~0 + replicate, y$samples)
# #### FULL MODEL ADJUST AND INTERESTING VARIABLES
# mod1 <- model.matrix(~0 + replicate + treatment, data=y$samples)
# 



#### svobj = sva(edata,mod,mod0,n.sv=n.sv)

svobj <- svaseq(svaob, mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv[c(1,2)])) #you can also specify to not use all sva, just 1,2, etc.


cleaned_count <- as.data.frame(cleanY(cpm(y[,c(1:6)]), mod1, svobj$sv[c(1,2)])) #you can also specify to not use all sva, just 1,2, etc.

cleaned_count1 <- as.data.frame(cleanY(norm_samples, mod1, svobj$sv[c(1,2)])) #you can also specify to not use all sva, just 1,2, etc.


cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.

cleaned_count <- as.data.frame(cleanY(as.matrix(y), mod1, svobj$sv))

cleaned_count <- as.data.frame(cleanY(norm_samples, mod1, svobj$sv)) 
cleaned_count <- as.data.frame(cleanY(norm_samples, mod1, svobj$sv[c(1,2)])) 

#exon

#combat and sva
adjustedc <- ComBat(cleaned_count, batch=sampleTable$replicate)
cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 

#sva and combat
adjustedc <- ComBat(counts, batch=sampleTable$replicate)
svobj <- svaseq(adjustedc, mod1, mod0) 
cleaned_count <- as.data.frame(cleanY(adjustedc, mod1, svobj$sv)) 

#sva alone
svaob<-as.matrix(counts)
svobj <- svaseq(svaob, mod1, mod0) 
cleaned_count <- as.data.frame(cleanY(svaob, mod1, svobj$sv)) 

#log_cleaned_count <- log2(cleaned_count)


### COMBAT 

adjustedc <- ComBat(y$counts, batch=y$samples$replicate)
#adjustedlcpm <- ComBat(lcpm_pre, batch=y$samples$replicate)
adjusted <- ComBat(as.matrix(cleaned_count), batch=y$samples$replicate)
adjusted <- ComBat(as.matrix(cleaned_count), batch=samples_vector_file$replicate)
adjusted1 <- ComBat(as.matrix(cleaned_count1), batch=samples_vector_file$replicate)
adjusted2 <- ComBat(as.matrix(adjusted), batch=y$samples$times)


#### ,MICROARRAY ANALYSIS FINAL DECISION SVA AND THEN COMBAT

svaob<-y$counts

#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod0 <- model.matrix(~1, y$samples)

#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + treatment, y$samples)

svobj <- svaseq(svaob, mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv[c(1,2)])) #you can also specify to not use all sva, just 1,2, etc.


#### COMBAT 
adjusted <- ComBat(as.matrix(cleaned_count), batch=y$samples$replicate)

####### RNA SEQ DECISION: SVA AND COMBAT ALL SV ON 8 AND 24 HOURS SEPARATELY

#### ,MICROARRAY ANALYSIS FINAL DECISION SVA AND THEN COMBAT 24 HOURS ALL
svaob<-y$counts

#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod0 <- model.matrix(~1, y$samples)

#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + treatment, y$samples)

svobj <- svaseq(svaob, mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


#### 24 HOURS
cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.


#### COMBAT 
adjusted <- ComBat(as.matrix(cleaned_count), batch=y$samples$replicate)

#### ,MICROARRAY ANALYSIS FINAL DECISION SVA AND THEN COMBAT 24 HOURS MKC VS DMSO
svaob<-as.matrix(y$counts[,c(1:3,7:9)])

#### NuLL MODEL WHEN WE DON'T KNOW THE BATCH
mod0 <- model.matrix(~1,y$samples[c(1:3,7:9),])

#### FULL MODEL ADJUST AND INTERESTING VARIABLES
mod1 <- model.matrix(~0 + treatment, y$samples[c(1:3,7:9),])

svobj <- svaseq(svaob, mod1, mod0) 

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}


#### 24 HOURS
cleaned_count <- as.data.frame(cleanY(cpm(y)[,c(1:3,7:9)], mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.
#cleaned_count <- as.data.frame(cleanY(cpm(y)[,c(1:3,7:9)], mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.


#### COMBAT 
adjusted <- ComBat(as.matrix(cleaned_count), batch=y$samples$replicate[c(1:3,7:9)])
