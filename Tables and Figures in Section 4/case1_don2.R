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
datalocation<- "Enter data location here"

#----------------------------------------------------------
datafilename1<- "4.D2UN_Uninfected CD4s.txt"
datafilename2<- "5.D2F4_HSA(hi).txt"
drop_list<- c('CD3','CD19','HSA')
pval.cutoff<-0.05
TRUH<-TRUE
Others<-!(TRUH)
#--------------------------------------------------------------

D_UN = read.csv(paste(datalocation,datafilename1,sep='/'),header=T,sep='\t')#change the name here
D_V = read.csv(paste(datalocation,datafilename2,sep='/'),header=T,sep='\t')#change the name here
t1<- D_UN
t2<- D_V
for (i in 1:(ncol(D_UN))){
  t1[which(D_UN[,i]<0),i]=0
  t2[which(D_V[,i]<0),i]=0
}
D_UN<-t1
D_V<- t2

index<-c(2:12,14:39,13)#putting HSA at the end
UN<-D_UN[,index];#change the name here
V<- D_V[,index];#change the name here

#2. Tests for Remodeling
UN<-D_UN[,2:39];#drop Event Type
V<- D_V[,2:39];#drop Event Type
index1<- which(names(UN) %in% drop_list)
index2<- which(names(V) %in% drop_list)
V.sub<-data.matrix(V[,-index2]);
UN.sub<-data.matrix(UN[,-index1]);
n.u<-dim(UN.sub)[1]
n.v<-dim(V.sub)[1]

if(Others){
  #---------------- Energy tests ------------------
  X <- rbind(UN.sub,V.sub)
  set.seed(1234)
  test1.energy<-eqdist.etest(X, sizes=c(n.u, n.v), distance=FALSE, method="original",R = 199)$p.value
  
  ## ---------------Crossmatch-------------------------
  set.seed(1)
  idx<-sample(1:n.u,ceiling(n.u/2),replace=FALSE)
  UN.sub.sampled<-UN.sub[idx,]
  X <- rbind(UN.sub.sampled,V.sub)
  n.u.new<-dim(UN.sub.sampled)[1]
  z <- c(rep(0,n.u.new),rep(1,n.v))
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
  test1.crossmatch<-crossmatchtest(z,dis)$approxpval
  rm(list=c('dis','out'))
  
  #-----------------gtests-------------------------------
  d <- dist(X)
  d <- as.matrix(d)
  G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
  rm(list=c('d','X'))
  mst<-minimum.spanning.tree(G)
  rm('G')
  E<-get.edgelist(mst,names=FALSE)
  out.gtests<-g.tests(E,1:n.u,(n.u+1):(n.u+n.v),test.type="all")
  test1.ec<-out.gtests$original$pval.approx
  test1.gec<-out.gtests$generalized$pval.approx
  test1.wec<-out.gtests$weighted$pval.approx
  test1.mtec<-out.gtests$maxtype$pval.approx
  gc()
  
  save(list=c('test1.energy','test1.crossmatch','out.gtests'),file="case1_don2_others.RData")
  
}

#----------------- Proposed ---------------------------
if(TRUH){
  out<-truh(V.sub,UN.sub,B=500,fc=1.1,a=pval.cutoff)
  
  save('out',file="case1_don2_truh.RData")
}

