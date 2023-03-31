library(mixOmics)

X<-t(readRDS("/home/bianca/Downloads/OCT_22_DIF_OMICS/rna_bianca/clean_DeRna_8"))

Y<-t(readRDS("/home/bianca/Downloads/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/clean_de_prot8"))

X<-t(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_DeRna_8"))
X<-X[rownames(Y),]
Y<-t(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/clean_de_prot8"))

X1<-t(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/rna_bianca/clean_DeRna_8"))
X<-X[rownames(Y),]Y<-t(readRDS("/home/ben/Desktop/OCT_22_DIF_OMICS/proteomics_cleaned_bianca/clean_de_prot8"))



samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
rownames(samples_info)=rownames(X)
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
colnames(samples_info)=c("condition")
rownames(samples_info)=rownames(X)
samples_info$replicate<-c(rep(c(1,2,3),2))

Y2<-Y*2 +5
spls.res <- spls(X = X, Y = Y2, ncomp = 5, mode = 'regression')
spls.res <- spls(X = X, Y = Y2, ncomp = 2, mode = 'regression')
spls.res <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
spls.res <- spls(X = X, Y = Y2,  mode = 'canonical')

# repeated CV tuning of component count
perf.spls.res <- perf(spls.res, validation = 'Mfold',
                        folds = 5, nrepeat = 5) 

plot(perf.spls.res, criterion = 'Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <-c(seq(100,500,4))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(seq(20,100,3))

# set range of test values for number of variables to use from X dataframe
list.keepX <-c(seq(20,100,4))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(seq(20,100,4))

# set range of test values for number of variables to use from X dataframe
list.keepX <-c(seq(10,50,4))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(seq(10,50,4))

### mode can also be regression if X explains Y
tune.spls.res <- tune.spls(X, Y, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 5, 
                             mode = 'canonical', measure = 'cor') 

tune.spls.res <- tune.spls(X, Y, ncomp = 2,
                           test.keepX = list.keepX,
                           test.keepY = list.keepY,
                           nrepeat = 1, folds = 5,
                           mode = 'regression', measure = 'cor')
plot(tune.spls.res)

tune.spls.res$choice.keepX

tune.spls.res$choice.keepY



# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.res$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.res$choice.keepY

optimal.ncomp <-  length(optimal.keepX)

# use all tuned values from above
final.spls.res <- spls(X, Y, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "canonical") 
final.spls.res <- spls(X, Y, ncomp = optimal.ncomp, 
                       keepX = optimal.keepX,
                       keepY = optimal.keepY,
                       mode = "regression") 
# explanitory approach being used, hence use regression when mode ="regression"




XimportantL<-selectVar(final.spls.res, comp = 1)$X$value
YimportantL<-selectVar(final.spls.res, comp = 1)$X$value
Ximportant<-selectVar(final.spls.res, comp = 1)$X$name
Yimportant<-selectVar(final.spls.res, comp = 1)$Y$name

XYimportant<-c(Ximportant, Yimportant)
saveRDS(Ximportant,"/home/ben/Desktop/mrna_proteins_spls_on_de8_rnasonly")
saveRDS(Yimportant,"/home/ben/Desktop/mrna_proteins_spls_on_de8_proteinsonly")
saveRDS(XYimportant,"/home/ben/Desktop/mrna_proteins_spls_on_de8")

vip.final.spls.res <- vip(final.spls.res)
X1<-vip.final.spls.res[!rowSums(vip.final.spls.res==0),]

head(vip.final.spls.res[selectVar(final.spls.res, comp = 1)$X$name,1])

perf.final.spls.res <- perf(final.spls.res, validation = 'Mfold', folds = 5, nrepeat = 5)
# Extract stability
stab.final.spls.res.comp1 <- perf.final.spls.res$features$stability.X$comp1
# Averaged stability of the X selected features across CV runs, as shown in Table
stab.final.spls.res.comp1[1:optimal.keepX[1]]

# We extract the stability measures of only the variables selected in final.spls.res
extr.stab.final.spls.res.Xcomp1 <- stab.final.spls.res.comp1[selectVar(final.spls.res, 
                                                                comp =1)$X$name]
extr.stab.final.spls.res.Xcomp2 <- stab.final.spls.res.comp1[selectVar(final.spls.res, 
                                                                      comp =2)$X$name]

extr.stab.final.spls.res.Ycomp1 <- stab.final.spls.res.comp1[selectVar(final.spls.res, 
                                                                      comp =1)$Y$name]
extr.stab.final.spls.res.Ycomp2 <- stab.final.spls.res.comp1[selectVar(final.spls.res, 
                                                                      comp =2)$Y$name]

plotIndiv(final.spls.res, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = samples_info$condition, # colour by time group
          pch = as.factor(samples_info$replicate), 
          col.per.group = color.mixo(1:2), 
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Replicate')

plotIndiv(final.spls.res, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = samples_info$condition, # colour by time group
          pch = as.factor(samples_info$replicate), 
          col.per.group = color.mixo(1:2), 
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Replicate')

plotIndiv(final.spls.res, ind.names = FALSE,
          rep.space = "XY-variate", # plot in Y-variate subspace
          group = samples_info$condition, # colour by time group
          pch = as.factor(samples_info$replicate), 
          col.per.group = color.mixo(1:2), 
          legend = TRUE, legend.title = 'Treatment', legend.title.pch = 'Replicate', title="mRNA and proteins")

## if you have more than 3 components
col.tox <- color.mixo(as.numeric(as.factor(samples_info$condition))) # create set of colours
plotIndiv(final.spls.res, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          axes.box = "both", col = col.tox, style = '3d')

### arrow plot
plotArrow(final.spls.res, ind.names = FALSE,
          group = samples_info$condition, # colour by time group
          #pch = as.factor(samples_info$replicate),
          col.per.group = color.mixo(1:2),
          legend.title = 'Treatment')


# form new perf() object which utilises the final model
perf.spls.res2 <- perf(final.spls.res, 
                        folds = 5, nrepeat = 10, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = FALSE)

# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(perf.spls.res2$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.res2$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))


plotVar(final.spls.res, cex = c(3,4), var.names = c(FALSE, TRUE))
plotVar(final.spls.res, cex = c(3,4), var.names = c(TRUE, TRUE))


color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

X11() # To open a new window for Rstudio
network(final.spls.res, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS res Toxicity Case Study Network Plot')


cim(final.spls.res, comp = 1:2, xlab = "proteins", ylab = "Genes")
