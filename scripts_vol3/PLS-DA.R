library(mixOmics)

dep<-readRDS("/home/bianca/Desktop/de_p_8_normal_values")
der<-readRDS("/home/bianca/Desktop/de_r_8_normal_values")
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
rownames(samples_info)=colnames(dep)
samples_info=as.data.frame(c(rep("DMSO",3), rep("MKC",3)))
colnames(samples_info)=c("condition")
rownames(samples_info)=colnames(dep)

########################### PCA - DA
# Biological question
# I am analysing a single data set (e.g. transcriptomics data) 
# and I would like to classify my samples into known groups 
# and predict the class of new samples. In addition, I am interested 
# in identifying the key variables that drive such discrimination.

X <- dep
Y <- samples_info$condition 
summary(Y)

dim(X); length(Y)

MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) # 1 Run the method
plotIndiv(MyResult.splsda)                          # 2 Plot the samples
plotVar(MyResult.splsda)                            # 3 Plot the variables
selectVar(MyResult.splsda, comp=1)$name             

MyResult.plsda <- plsda(X,Y) # 1 Run the method
plotIndiv(MyResult.plsda)    # 2 Plot the samples
plotVar(MyResult.plsda)      # 3 Plot the variables

plotIndiv(MyResult.splsda, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, star = TRUE, title = 'sPLS-DA on SRBCT',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

plotVar(MyResult.splsda, var.names=FALSE)
plotVar(MyResult.plsda, cutoff=0.7)

background <- background.predict(MyResult.splsda, comp.predicted=2,
                                 dist = "max.dist") 
plotIndiv(MyResult.splsda, comp = 1:2, group = srbct$class,
          ind.names = FALSE, title = "Maximum distance",
          legend = TRUE,  background = background)

auc.plsda <- auroc(MyResult.splsda)


MyResult.splsda2 <- splsda(X,Y, ncomp=3, keepX=c(15,10,5))
selectVar(MyResult.splsda2, comp=1)$value

plotLoadings(MyResult.splsda2, contrib = 'max', method = 'mean')
plotIndiv(MyResult.splsda2, style="3d")

MyResult.plsda2 <- plsda(X,Y, ncomp=10)
set.seed(30) # for reproducibility in this vignette, otherwise increase nrepeat
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 3, 
                     progressBar = FALSE, nrepeat = 10) # we suggest nrepeat = 50

# type attributes(MyPerf.plsda) to see the different outputs
# slight bug in our output function currently see the quick fix below
#plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# quick fix
matplot(MyPerf.plsda$error.rate$BER, type = 'l', lty = 1, 
        col = color.mixo(1:3), 
        main = 'Balanced Error rate')
legend('topright', 
       c('max.dist', 'centroids.dist', 'mahalanobis.dist'), 
       lty = 1,
       col = color.mixo(5:7))


MyPerf.plsda

list.keepX <- c(5:10,  seq(20, 100, 10))
list.keepX # to output the grid of values tested

tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 3, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 3, dist = 'max.dist', progressBar = FALSE,
                                 measure = "BER", test.keepX = list.keepX,
                                 nrepeat = 10)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
select.keepX

plot(tune.splsda.srbct, col = color.jet(ncomp))

MyResult.splsda.final <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX)
plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="SPLS-DA, Final result")
