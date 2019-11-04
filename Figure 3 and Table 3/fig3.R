
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
library(latex2exp)

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
experiment.d2.1<-function(n.u,n.v,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  
  result<-foreach(r = 1:reps,.packages=c("MASS","crossmatch","gTests","energy","igraph"),.export="nearest")%dopar%{
    
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
    test1.energy<-1*(eqdist.etest(X, sizes=c(n.u, n.v), distance=FALSE, method="original",R = 199)$p.value<pval.cutoff)
    
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
    test1.crossmatch<-1*(crossmatchtest(z,dis)$approxpval<pval.cutoff)
    
    #-----------------gtests-------------------------------
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test1.ec<-1*(out.gtests$original$pval.approx<pval.cutoff)
    test1.gec<-1*(out.gtests$generalized$pval.approx<pval.cutoff)
    test1.wec<-1*(out.gtests$weighted$pval.approx<pval.cutoff)
    test1.mtec<-1*(out.gtests$maxtype$pval.approx<pval.cutoff)
    
    #----------------- Proposed ---------------------------
    out.nearest<- matrix(0,n.v,3)
    for(i in 1:dim(V)[1]){
      out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
    }
    test1.nearest<-mean(out.nearest[,3])
    return(list("test1.energy"=test1.energy, "test1.crossmatch"=test1.crossmatch,"test1.ec"=test1.ec,
                "test1.gec"=test1.gec,"test1.wec"=test1.wec,"test1.mtec"=test1.mtec,
                "test1.nearest"=test1.nearest,"U"=U,"V"=V))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  out<-cbind(mean(sapply(1:reps,function(i) result[[i]]$test1.energy)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.crossmatch)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.ec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.gec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.wec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.mtec)))
  test1.nearest<- sapply(1:reps,function(i) result[[i]]$test1.nearest)
  return(list("tests"=out,"nearest"=test1.nearest,"U"=result[[1]]$U,
              "V"=result[[1]]$V))
}
#----------------------------------------------------------
experiment.d2.2<-function(n.u,n.v,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  
  result<-foreach(r = 1:reps,.packages=c("MASS","crossmatch","gTests","energy","igraph"),.export="nearest")%dopar%{
    
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
    test1.energy<-1*(eqdist.etest(X, sizes=c(n.u, n.v), distance=FALSE, method="original",R = 199)$p.value<pval.cutoff)
    
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
    test1.crossmatch<-1*(crossmatchtest(z,dis)$approxpval<pval.cutoff)
    
    #-----------------gtests-------------------------------
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test1.ec<-1*(out.gtests$original$pval.approx<pval.cutoff)
    test1.gec<-1*(out.gtests$generalized$pval.approx<pval.cutoff)
    test1.wec<-1*(out.gtests$weighted$pval.approx<pval.cutoff)
    test1.mtec<-1*(out.gtests$maxtype$pval.approx<pval.cutoff)
    
    #----------------- Proposed ---------------------------
    out.nearest<- matrix(0,n.v,3)
    for(i in 1:dim(V)[1]){
      out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
    }
    test1.nearest<-mean(out.nearest[,3])
    return(list("test1.energy"=test1.energy, "test1.crossmatch"=test1.crossmatch,"test1.ec"=test1.ec,
                "test1.gec"=test1.gec,"test1.wec"=test1.wec,"test1.mtec"=test1.mtec,
                "test1.nearest"=test1.nearest,"U"=U,"V"=V))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  out<-cbind(mean(sapply(1:reps,function(i) result[[i]]$test1.energy)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.crossmatch)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.ec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.gec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.wec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.mtec)))
  test1.nearest<- sapply(1:reps,function(i) result[[i]]$test1.nearest)
  return(list("tests"=out,"nearest"=test1.nearest,"U"=result[[1]]$U,
              "V"=result[[1]]$V))
}
#----------------------------------------------------------
experiment.d2.3<-function(n.u,n.v,pval.cutoff,reps){
  
  d.u<-2
  d.v<-2
  Sigma.u = diag(d.u)
  Sigma.v = diag(d.v)
  
  cl <- makeCluster(8) 
  registerDoParallel(cl)
  
  result<-foreach(r = 1:reps,.packages=c("MASS","crossmatch","gTests","energy","igraph"),.export="nearest")%dopar%{
    
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
    test1.energy<-1*(eqdist.etest(X, sizes=c(n.u, n.v), distance=FALSE, method="original",R = 199)$p.value<pval.cutoff)
    
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
    test1.crossmatch<-1*(crossmatchtest(z,dis)$approxpval<pval.cutoff)
    
    #-----------------gtests-------------------------------
    d <- as.matrix(d)
    G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
    mst<-minimum.spanning.tree(G)
    E<-get.edgelist(mst,names=FALSE)
    out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
    test1.ec<-1*(out.gtests$original$pval.approx<pval.cutoff)
    test1.gec<-1*(out.gtests$generalized$pval.approx<pval.cutoff)
    test1.wec<-1*(out.gtests$weighted$pval.approx<pval.cutoff)
    test1.mtec<-1*(out.gtests$maxtype$pval.approx<pval.cutoff)
    
    #----------------- Proposed ---------------------------
    out.nearest<- matrix(0,n.v,3)
    for(i in 1:dim(V)[1]){
      out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
    }
    test1.nearest.usd<-mean(out.nearest[,3])
    test1.nearest.sd<-sqrt(n.v)*mean(out.nearest[,3])/sd(out.nearest[,3])
    sd1<-sd(out.nearest[,1])
    sd2<-sd(out.nearest[,2])
    temp1<-(sd1/sd2)^2
    temp2<-(sd2/sd1)^2
    w1<-temp1/(temp1+temp2)
    w2<-1-w1
    test1new.nearest.usd<-mean(w1*out.nearest[,1]-w2*out.nearest[,2])
    test1new.nearest.sd<-sqrt(n.v)*mean(w1*out.nearest[,1]-w2*out.nearest[,2])/sd(w1*out.nearest[,1]-w2*out.nearest[,2])
    
    return(list("test1.energy"=test1.energy, "test1.crossmatch"=test1.crossmatch,"test1.ec"=test1.ec,
                "test1.gec"=test1.gec,"test1.wec"=test1.wec,"test1.mtec"=test1.mtec,
                "test1.nearest.usd"=test1.nearest.usd,
                "test1.nearest.sd"=test1.nearest.sd,
                "test1new.nearest.usd"=test1new.nearest.usd,
                "test1new.nearest.sd"=test1new.nearest.sd,
                "U"=U,"V"=V))
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  out<-cbind(mean(sapply(1:reps,function(i) result[[i]]$test1.energy)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.crossmatch)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.ec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.gec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.wec)),
             mean(sapply(1:reps,function(i) result[[i]]$test1.mtec)))
  test1.nearest.usd<- sapply(1:reps,function(i) result[[i]]$test1.nearest.usd)
  test1.nearest.sd<- sapply(1:reps,function(i) result[[i]]$test1.nearest.sd)
  test1new.nearest.usd<- sapply(1:reps,function(i) result[[i]]$test1new.nearest.usd)
  test1new.nearest.sd<- sapply(1:reps,function(i) result[[i]]$test1new.nearest.sd)
  return(list("tests"=out,"nearest.usd"=test1.nearest.usd,
              "nearest.sd"=test1.nearest.sd,"nearestnew.usd"=test1new.nearest.usd,
              "nearestnew.sd"=test1new.nearest.sd,
              "U"=result[[1]]$U,
              "V"=result[[1]]$V))
}
#----------------------------------------------------------
#***********************************************************************************************
n.u<-2000
n.v<-500
reps<-100
out.d2.1<-experiment.d2.1(n.u,n.v,0.05,reps)

plotdata.1<-as.data.frame(rbind(out.d2.1$U,out.d2.1$V))
plotdata.1$Case<- as.factor(c(rep('Uinfected',n.u),rep('Infected',n.v)))

g1<-ggplot(plotdata.1, aes(x = V1, y = V2,color=Case,shape=Case)) +geom_point()+
  xlab(TeX("Marker 1 ($X_1$)"))+ylab(TeX("Marker 2 ($X_2$)"))+theme_bw()

#--------------------------------------------------------------------
out.d2.2<-experiment.d2.2(n.u,n.v,0.05,reps)

plotdata.2<-as.data.frame(rbind(out.d2.2$U,out.d2.2$V))
plotdata.2$Case<- as.factor(c(rep('Uinfected',n.u),rep('Infected',n.v)))

g2<-ggplot(plotdata.2, aes(x = V1, y = V2,color=Case,shape=Case)) +geom_point()+
  xlab(TeX("Marker 1 ($X_1$)"))+ylab(TeX("Marker 2 ($X_2$)"))+theme_bw()

#--------------------------------------------------------------------
out.d2.3<-experiment.d2.3(n.u,n.v,0.05,reps)

plotdata.3<-as.data.frame(rbind(out.d2.3$U,out.d2.3$V))
plotdata.3$Case<- as.factor(c(rep('Uinfected',n.u),rep('Infected',n.v)))

g3<-ggplot(plotdata.3, aes(x = V1, y = V2,color=Case,shape=Case)) +geom_point()+
  xlab(TeX("Marker 1 ($X_1$)"))+ylab(TeX("Marker 2 ($X_2$)"))+theme_bw()

#--------------------------------------------------------------------

ggarrange(g1, g2, g3, ncol=3, nrow=1, common.legend = TRUE, legend="top")
