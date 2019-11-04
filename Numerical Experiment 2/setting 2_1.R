
library(MASS)
library(crossmatch)
library(gTests)
library(energy)
library(ggplot2)
library(gridExtra)
library(igraph)
library(foreach)
library(doParallel)
library(iterators)
library(clusterGeneration)
library(lcmix)
library(CVTuningCov)

source('lib.R')

reps=100
pval.cutoff=0.05
BB=200

#----------------------------------------------------------
experiment.1<-function(n.u,n.v,d,pval.cutoff,reps,nB){
  
  Sigma.1 = AR1(d,0.7)
  Sigma.2 = AR1(d,-0.9)
 
  test.truh<-matrix(0,reps,2)
  test.others<-matrix(0,reps,6)
  
  for(r in 1:reps){
    set.seed(r)
    p<-runif(n.u,0,1)
    set.seed(10*r)
    U<- (p<=0.5)*rmvgamma(n.u,shape=runif(d,5,5),rate=rep(1,d),Sigma.1)+(p>0.5)*rmvexp(n.u,rate=rep(1,d),Sigma.2)
    set.seed(50*r)
    V<-rmvexp(n.v,rate=rep(1,d),Sigma.2)
    #---------------- Energy tests ------------------
    X <- rbind(U,V)
    set.seed(1234)
    test.others[r,1]<-eqdist.etest(X, sizes=c(n.u, n.v), distance=FALSE, method="original",R = 199)$p.value

    ## ---------------Crossmatch-------------------------

    z <- c(rep(0,n.u),rep(1,n.v))
    #Rank based Mahalanobis distance between each pair:
    X <- as.matrix(X)
    n <- dim(X)[1]
    k <- dim(X)[2]
    for (j in 1:k) X[,j] <- rank(X[,j])
    cv <- cov(X)
    vuntied <- var(1:n)
    rat <- sqrt(vuntied/diag(cv))
    cv <- diag(rat)%*%cv%*%diag(rat)
    out <- matrix(NA,n,n)
    icov <- ginv(cv)
    for (i in 1:n) out[i,] <- mahalanobis(X,X[i,],icov,inverted=TRUE)
    dis <- out
    test.others[r,2]<-crossmatchtest(z,dis)$approxpval
    rm(list=c('dis','out'))

    #-----------------gtests-------------------------------
    dd <- dist(X)
    dd <- as.matrix(dd)
    G<-graph.adjacency(dd,mode="undirected",weighted=TRUE)
    rm(list=c('dd','X'))
    mst<-minimum.spanning.tree(G)
    rm('G')
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test.others[r,3]<-out.gtests$original$pval.approx
    test.others[r,4]<-out.gtests$generalized$pval.approx
    test.others[r,5]<-out.gtests$weighted$pval.approx
    test.others[r,6]<-out.gtests$maxtype$pval.approx
    gc()
    #----------------- Proposed ---------------------------
    out<-truh(V,U,nB,fc=1,a=pval.cutoff)
    test.truh[r,]<-c(out$teststat,out$pval)
    
  }
  return(list("truh"=test.truh, "test.others"=test.others))
}
#----------------------------------------------------------

#***********************************************************************************************
d<- 5 #rho = 0.1
out.d1.1<-experiment.1(500,50,d,pval.cutoff,reps,BB)
out.d1.2<-experiment.1(2000,200,d,pval.cutoff,reps,BB)

#***********************************************************************************************
d<- 15 # rho = 0.1
out.d2.1<-experiment.1(500,50,d,pval.cutoff,reps,BB)
out.d2.2<-experiment.1(2000,200,d,pval.cutoff,reps,BB)

#***********************************************************************************************
d<- 30 # rho = 0.1
out.d3.1<-experiment.1(500,50,d,pval.cutoff,reps,BB)
out.d3.2<-experiment.1(2000,200,d,pval.cutoff,reps,BB)

rejrate.1<-matrix(0,1,3)
rejrate.2<-matrix(0,1,3)

for(i in 1:6){

  rejrate.1[i,]<-c(mean(1*(out.d1.1$test.others[,i])<=pval.cutoff),
                      mean(1*(out.d2.1$test.others[,i])<=pval.cutoff),
                      mean(1*(out.d3.1$test.others[,i])<=pval.cutoff))

  rejrate.2[i,]<-c(mean(1*(out.d1.2$test.others[,i])<=pval.cutoff),
                      mean(1*(out.d2.2$test.others[,i])<=pval.cutoff),
                      mean(1*(out.d3.2$test.others[,i])<=pval.cutoff))
}


save.image(paste(getwd(),'/setting2_1.RData',sep=""))


