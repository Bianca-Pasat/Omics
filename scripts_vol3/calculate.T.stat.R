#Two group data (not paired): calculate the wilcoxon rank sum (WRS) statistic
#data are the expression matrix:n_iso*n_s
#tag are the permutation matrix:n_perm*n_s
#return a matrix: n_iso*n_perm
calculate.WRS.stat=function(data, tag)
{
  data=apply(data,2,as.numeric) #this is important, since data is a dataframe object, need to convert to numeric
  n_iso=nrow(data)
  n_perm=nrow(tag)
  n_s=ncol(data)
  n1=sum(tag[1,]==1)
  
  rank_sum=matrix(0,nrow=n_iso,ncol=n_perm) #this matrix will take the rank sum
  
  for(i in 1:n_perm)
  {
    #print(tag[i,])
    mat=rbind(data,tag[i,])
    #print(mat[(n_iso+1),]==1)
    #print(dim(mat))
    y1=mat[1:n_iso,mat[(n_iso+1),]==1]
    y2=mat[1:n_iso,mat[(n_iso+1),]==2]
    
    r=t(apply(cbind(y1,y2),1,rank))
    
    rank_sum[,i]=rowSums(r[,(1:n1)])
  }
  
  return( rank_sum-n1*(n_s+1)/2 )
}

#Two paired group data: calculate the wilcoxon signed rank(WRS) statistic
#data are the expression matrix:n_iso*n_s
#tag are the permutation matrix:n_perm*n_s
#return a matrix: n_iso*n_perm
calculate.WSR.stat=function(data, tag)
{
  data=apply(data,2,as.numeric) #this is important, since data is a dataframe object, need to convert to numeric
  n_iso=nrow(data)  
  n_perm=nrow(tag)
  n_s=ncol(data)
  
  WSR=matrix(0,nrow=n_iso,ncol=n_perm) #this matrix will take the WSR stat  
  
  for(i in 1:n_perm)
  {
    mat=rbind(data,tag[i,])
    
    y1=mat[1:n_iso,mat[(n_iso+1),]==1]
    y2=mat[1:n_iso,mat[(n_iso+1),]==2]
    
    x=y1-y2 #the difference
    
    for(j in 1:n_iso)
    {
      xj=x[j,]
      ZEROES=any(xj==0)	  
      if(ZEROES) xj=xj[xj != 0]
      n_a=length(xj)
      r=rank(abs(xj))
      WSR[j,i]=sum(r[xj>0])-(n_a+1)*n_a/4
    } 
  }
  
  return( WSR )
}

#Multi-group data: calculate the Kruskal-Wallis(KW) statistic
#data is the expression matrix:n_iso*n_s
#tag is the group tag vector
#perm_index is the index for the permutation matrix:n_perm*n_s
#return a matrix: n_iso*n_perm
calculate.KW.stat=function(data,tag,perm_index)
{
  data=apply(data,2,as.numeric) #this is important, since data is a dataframe object, need to convert to numeric
  n_iso=nrow(data)  
  n_perm=nrow(perm_index)
  n_s=ncol(data)
  
  KW=matrix(0,nrow=n_iso,ncol=n_perm) #this matrix will take the WSR stat  
  
  for(i in 1:n_iso)
  {  
    r = rank(data[i,])
    TIES = table(data[i,])
    
    for(j in 1:n_perm)
    {
      g = factor(tag[perm_index[j,]])
      STATISTIC = sum(tapply(r, g, "sum")^2/tapply(r, g, "length"))
      KW[i,j] = ((12 * STATISTIC/(n_s * (n_s + 1)) - 3 * (n_s + 1))/(1 - sum(TIES^3 - TIES)/(n_s^3 - n_s)))
    }
  }
  
  return(KW)
}

#Quantitative outcome: calculate the Spearman's correlation coefficient (SCC) statistic
#data is the expression matrix:n_iso*n_s
#outcome is the quantitative outcome vector
#perm_index is the index for the permutation matrix:n_perm*n_s
#return a matrix: n_iso*n_perm
calculate.SCC.stat=function(data,outcome,perm_index)
{
  data=t(apply(data,2,as.numeric))
  #this is important, since data is a dataframe object, need to convert to numeric
  #also need to take transpose to calculate the SCC
  
  n_perm=dim(perm_index)[1]
  n_s=dim(perm_index)[2]
  n_iso=dim(data)[2] #data matrix is transposed
  
  outcome=matrix(outcome[t(perm_index)],nrow=n_perm,byrow=T)
  
  cor_matrix=matrix(0, nrow=n_iso, ncol=n_perm)
  
  for(i in 1:n_perm)
  {
    cor_matrix[,i]=cor(data, outcome[i,], method='spearman')
  }
  
  return(cor_matrix)
}

#Survival outcome: calculate the Cox proportional hazard model score statistic (CS) 
#data is the expression matrix:n_iso*n_s
#outcome is the quantitative outcome vector
#perm_index is the index for the permutation matrix:n_perm*n_s
#return a matrix: n_iso*n_perm
calculate.CS.stat <- function(data, time, censor ,perm_index)
  #mu is the fitted values (can be a matrix); y is the survival time vector; gamma is the censored status vector 
{
  data=apply(data,2,as.numeric) #this is important, since data is a dataframe object, need to convert to numeric
  n_iso=nrow(data)  
  n_perm=nrow(perm_index)
  n_s=ncol(data)
  
  mu=apply(data,1,rank,ties.method="random")
  
  score=matrix(0, nrow=n_iso)
  for(i in 1:n_perm)
  {
    y=time[perm_index[i,]]
    gamma=censor[perm_index[i,]]
    
    # find the index matrix
    Dn <- sum(gamma == 1)
    Dset <- c(1 : ncol(mu))[gamma == 1]		# the set of observed
    ind <- matrix(0, ncol(mu), Dn)
    
    # get the matrix
    for (i in 1 : Dn)
    {
      ind[y > y[Dset[i]] - 1e-8, i] <- 1 / sum(y > y[Dset[i]] - 1e-8)
    }
    ind.sums <- rowSums(ind)
    x.ind <- mu %*% ind
    
    # get the derivatives
    dev1 <- mu %*% (gamma - ind.sums)
    dev2 <- (mu * mu) %*% ind.sums - rowSums(x.ind * x.ind)
    
    score[,i] <- dev1 / (sqrt(-dev2) + SMALL.VAL)
  }
  
  return(score)
}
