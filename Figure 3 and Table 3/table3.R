
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
library(ggpubr)

source('lib.R')
  
#----------------------------------------------------------
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x +
                 theme(legend.position = "none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl), 
                                            legend,ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend, ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}
#----------------------------------------------------------
experiment.d2.1<-function(n.u,n.v,nB,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  test.truh<-matrix(0,reps,1)
  test.others<-matrix(0,reps,6)
  
 for(r in 1:reps){
    
    set.seed(r)
    p<-runif(n.u,0,1)
    set.seed(10*r)
    U<- (p<=0.3)*mvrnorm(n.u,c(0,0),Sigma.u)+(p>0.3 & p<=0.6)*mvrnorm(n.u,c(0,-4),Sigma.u)+(p>0.6)*mvrnorm(n.u,c(4,-2),Sigma.u)
    set.seed(20*r)
    q<-runif(n.v,0,1)
    set.seed(50*r)
    V<-(q<=0.3)*mvrnorm(n.v,c(0,0),Sigma.v)+(q>0.3 & q<=0.6)*mvrnorm(n.v,c(0,-4),Sigma.v)+(q>0.6)*mvrnorm(n.v,c(4,-2),Sigma.v)
    
    #---------------- Energy tests ------------------
    X <- rbind(U,V)
    d <- dist(X)
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
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test.others[r,3]<-out.gtests$original$pval.approx
    test.others[r,4]<-out.gtests$generalized$pval.approx
    test.others[r,5]<-out.gtests$weighted$pval.approx
    test.others[r,6]<-out.gtests$maxtype$pval.approx
    
    #----------------- Proposed ---------------------------
    out<-truh(V,U,nB,fc=1,a=pval.cutoff)
    test.truh[r,1]<-out$pval
    
  }
  
  rejrate.mat<-c(colMeans(1*(test.others<pval.cutoff)),colMeans(1*(test.truh<0.05)))
  return(rejrate.mat)
}
#----------------------------------------------------------
experiment.d2.2<-function(n.u,n.v,nB,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  test.truh<-matrix(0,reps,1)
  test.others<-matrix(0,reps,6)
  
  for(r in 1:reps){
    
    set.seed(r)
    p<-runif(n.u,0,1)
    set.seed(10*r)
    U<- (p<=0.3)*mvrnorm(n.u,c(0,0),Sigma.u)+(p>0.3 & p<=0.6)*mvrnorm(n.u,c(0,-4),Sigma.u)+(p>0.6)*mvrnorm(n.u,c(4,-2),Sigma.u)
    set.seed(20*r)
    q<-runif(n.v,0,1)
    set.seed(50*r)
    V<-(q<=0.5)*mvrnorm(n.v,c(2,-2),Sigma.v)+(q>0.5)*mvrnorm(n.v,c(3,3),Sigma.v)
    
    #---------------- Energy tests ------------------
    X <- rbind(U,V)
    d <- dist(X)
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
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test.others[r,3]<-out.gtests$original$pval.approx
    test.others[r,4]<-out.gtests$generalized$pval.approx
    test.others[r,5]<-out.gtests$weighted$pval.approx
    test.others[r,6]<-out.gtests$maxtype$pval.approx
    
    #----------------- Proposed ---------------------------
    out<-truh(V,U,nB,fc=1,a=pval.cutoff)
    test.truh[r,1]<-out$pval
    
  }
  
  rejrate.mat<-c(colMeans(1*(test.others<pval.cutoff)),colMeans(1*(test.truh<0.05)))
  return(rejrate.mat)
}
#----------------------------------------------------------
experiment.d2.3<-function(n.u,n.v,nB,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  test.truh<-matrix(0,reps,1)
  test.others<-matrix(0,reps,6)
  
  for(r in 1:reps){
    set.seed(r)
    p<-runif(n.u,0,1)
    set.seed(10*r)
    U<- (p<=0.3)*mvrnorm(n.u,c(0,0),Sigma.u)+(p>0.3 & p<=0.6)*mvrnorm(n.u,c(0,-4),Sigma.u)+(p>0.6)*mvrnorm(n.u,c(4,-2),Sigma.u)
    set.seed(20*r)
    q<-runif(n.v,0,1)
    set.seed(50*r)
    V<-(q<=0.8)*mvrnorm(n.v,c(0,0),Sigma.v)+(q>0.8 & q<=0.9)*mvrnorm(n.v,c(0,-4),Sigma.v)+(q>0.9)*mvrnorm(n.v,c(4,-2),Sigma.v)
    
    #---------------- Energy tests ------------------
    X <- rbind(U,V)
    d <- dist(X)
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
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test.others[r,3]<-out.gtests$original$pval.approx
    test.others[r,4]<-out.gtests$generalized$pval.approx
    test.others[r,5]<-out.gtests$weighted$pval.approx
    test.others[r,6]<-out.gtests$maxtype$pval.approx
    
    #----------------- Proposed ---------------------------
    out<-truh(V,U,nB,fc=1,a=pval.cutoff)
    test.truh[r,1]<-out$pval
  }
  
rejrate.mat<-c(colMeans(1*(test.others<pval.cutoff)),colMeans(1*(test.truh<0.05)))
return(rejrate.mat)
}
#----------------------------------------------------------
#***********************************************************************************************
n.u<-2000
n.v<-500
reps<-100
nB<-200
out.d2.1<-experiment.d2.1(n.u,n.v,nB,0.05,reps)
#--------------------------------------------------------------------
out.d2.2<-experiment.d2.2(n.u,n.v,nB,0.05,reps)
#--------------------------------------------------------------------
out.d2.3<-experiment.d2.3(n.u,n.v,nB,0.05,reps)

