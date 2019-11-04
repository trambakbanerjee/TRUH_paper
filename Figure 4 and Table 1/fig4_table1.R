library(gTests)
library(igraph)


nearest<-function(y,U,n=100,d=10){
  b1=matrix(y,ncol=d,nrow=n,byrow=T)
  temp=rowsums(abs(b1-U))
  ind1=which.min(temp);
  y.u=U[ind1,];
  d1=temp[ind1];
  
  b2=matrix(y.u,ncol=d,nrow=n-1,byrow=T)
  U1=as.matrix(U[-c(ind1),]);
  temp=rowsums(abs(b2-U1))
  ind2=which.min(temp);
  u.u=U1[ind2,];
  d2=temp[ind2];
  return(c(d1,d2,d1-d2))
}

truh<-function(V,U,B,fc=1,a=0.05){
  
  n.u<-dim(U)[1]
  n.v<-dim(V)[1]
  d.u<-dim(U)[2]
  d.v<-dim(V)[2]
  
  out.nearest<- matrix(0,n.v,3)
  for(i in 1:dim(V)[1]){
    out.nearest[i,]<-nearest(V[i,],U,n.u,d.u)
  }
  fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
  teststat.truh<-fac*mean(out.nearest[,1]-1*out.nearest[,2])
  
  #------------ cutoff computation --------------------
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  
  X<-U
  set.seed(1)
  k.hat<-3
  clust.X <-clara(X,k=k.hat,metric="euclidean",samples=50,sampsize=500)$clustering
  mu<-table(clust.X)
  temp1<-diag(k.hat)
  w.v.b<-temp1[rep(1:k.hat,each=B),]
  s.v<- ceiling(w.v.b*n.v)
  Bp<-dim(s.v)[1]
  if(k.hat==1){
    Bp = B
  }
  
  result<-foreach(b = 1:Bp,.packages="Rfast",.export="nearest")%dopar%{
    
    list.u.b<-list()
    list.v.b<-list()
    if(k.hat>1){
      for(k in 1:k.hat){
        clust.idx<-which(clust.X==k)
        
        set.seed(k*b)
        idx1.v<-sample(clust.idx,s.v[b,k],replace=(s.v[b,k]>mu[k]))
        idx1.u<-clust.idx[!(clust.idx %in% idx1.v)]
        list.v.b[[k]]<- X[idx1.v,]
        list.u.b[[k]]<-X[idx1.u,]
      }
    } else{
      clust.idx<-which(clust.X==1)
      
      set.seed(b)
      idx1.v<-sample(clust.idx,s.v[b],replace=(s.v[b]>mu[1]))
      idx1.u<-clust.idx[!(clust.idx %in% idx1.v)]
      list.v.b[[1]]<- X[idx1.v,]
      list.u.b[[1]]<-X[idx1.u,]
      
    }
    if(d.u>1){
      VV.b<- do.call(rbind,list.v.b)
      UU.b<- do.call(rbind,list.u.b)
    } else {
      VV.b<-as.matrix(unlist(list.v.b))
      UU.b<-as.matrix(unlist(list.u.b))
    }
    out<- matrix(0,dim(VV.b)[1],2)
    for(i in 1:dim(VV.b)[1]){
      out[i,]<-nearest(VV.b[i,],UU.b,dim(UU.b)[1],d.u)[1:2]
    }
    
    fac<-(n.v^{1/d.u})*(d.u>1)+1*(d.u==1)
    dist.null<-fac*mean(fc*out[,1]-1*out[,2])
  }
  dist.null<-do.call(rbind,result)
  #*******************************************************************
  c1<-quantile(dist.null,a/2)
  c2<-quantile(dist.null,1-a/2)
  pval<-sum(1*(dist.null>teststat.truh))/Bp
  stopCluster(cl)
  registerDoSEQ()
  return(list("teststat"=teststat.truh,"c1"=c1,"c2"=c2,"pval"=pval,"null.dist"=dist.null))
  
}


L=500;
N=1000;
M=50;
nB=50;
pval.cutoff<-0.05
fr.A=matrix(0,L,2)
fr.B<-fr.A
truh.A=fr.A
truh.B<-fr.B
temp=matrix(0,M,1)
for (i in 1:L)
{
  set.seed(i)
  u=rnorm(N)
  v=rnorm(M)
  mu=c(rep(0,0.33*N),rep(4,0.34*N),rep(8,0.33*N));
  mu1=c(rep(2,0.33*N),rep(6,0.34*N),rep(10,0.33*N));
  U=2.5*mu+u;
  V=2.5*mu1+v;

  X <- rbind(U,V)
  d <- dist(t(X))
  d <- as.matrix(d)
  G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
  mst<-minimum.spanning.tree(G)
  E<-get.edgelist(mst,names=FALSE)
  out.gtests<-g.tests(E,1:N,(N+1):(N+M),test.type="all")
  fr.A[i,]=c(out.gtests$original$test.statistic,out.gtests$original$pval.approx)
  out<-truh(as.matrix(V),as.matrix(U),nB,fc=1,a=pval.cutoff)
  truh.A[i,]<-c(out$teststat,out$pval)
  print(i)
}

truh1=truh.A[,1]
fr1=fr.A[,1]
U1=U
V1=V

#------------------CASE B---------------------------------

for (i in 1:L)
{
  set.seed(i)
  u=rnorm(N)
  v=rnorm(M)
  mu=c(rep(0,0.33*N),rep(4,0.34*N),rep(8,0.33*N));
  mu1=c(rep(0,0.5*N),rep(8,0.5*N));
  U=2.5*mu+u;
  V=2.5*mu1+v;
  X <- rbind(U,V)
  d <- dist(t(X))
  d <- as.matrix(d)
  G<-graph.adjacency(d,mode="undirected",weighted=TRUE)
  mst<-minimum.spanning.tree(G)
  E<-get.edgelist(mst,names=FALSE)
  out.gtests<-g.tests(E,1:N,(N+1):(N+M),test.type="all")
  fr.B[i,]=c(out.gtests$original$test.statistic,out.gtests$original$pval.approx)
  out<-truh(as.matrix(V),as.matrix(U),nB,fc=1,a=pval.cutoff)
  truh.B[i,]<-c(out$teststat,out$pval)
  print(i)
}

save.image("fig2_table1.Rdata")


#----------------- Plots ---------------------------

par(mfrow=c(2,2))
plot(density(U1),main="Case A: Remodeling",xlab="",
     ylim=c(0,0.1),lwd="2",col=rgb(0,0,1,0.5),cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1.5)
points(density(V1),main="",xlab="",type="l",col=rgb(0,0,1,0.5),lwd=2,lty=2)
box(lwd=2,col="blue")


plot(density(U),main="Case B: No Remodeling, Only Preferential Infection",
     xlab="",ylim=c(0,0.1),lwd="2",col=rgb(1,0,0,0.5),cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1.5)
points(density(V),main="",xlab="",type="l",col=rgb(1,0,0,0.5),lwd=2,lty=2)
box(lwd=2,col="red")

t=truh.B[,1]
t1=truh1

hist(c(t,t1),freq=T,border="white",main="Distribution of TRUH statistics in cases A (blue) and B (red)",
     xlab="values",cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1.5)
hist(t,col=rgb(1,0,0,0.5),freq=T,add=T,breaks=1)
hist(t1,col=rgb(0,0,1,0.5),freq=T,add=T)

q=quantile(t,1-c(0.01,0.05,0.1,0.2));
power=q;

for (i in 1:length(q)){
power[i]=mean(as.numeric(t1>q[i]))
}

cutoff.truh=q;
power.truh=power;

rejectionrate.truh.B<-c(mean(1*(truh.B[,2]<0.01)),mean(1*(truh.B[,2]<0.05)),
                    mean(1*(truh.B[,2]<0.1)),mean(1*(truh.B[,2]<0.2)))

#-----------------------
#-----------------------

t=fr.B[,1]
t1=fr1

hist(c(t,t1),freq=T,border="white",main="Distribution of Edgecount statistics in cases A (blue) and B (red)",
     xlab="values",cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1.5)
hist(t,col=rgb(1,0,0,0.5),freq=T,add=T)
hist(t1,col=rgb(0,0,1,0.5),freq=T,add=T)

q=quantile(t,1-c(0.01,0.05,0.1,0.2));
power=q;

for (i in 1:length(q)){
power[i]=mean(as.numeric(t1>q[i]))
}

cutoff.ec=q;
power.ec=power;

rejectionrate.ec.B<-c(mean(1*(fr.B[,2]<0.01)),mean(1*(fr.B[,2]<0.05)),
                        mean(1*(fr.B[,2]<0.1)),mean(1*(fr.B[,2]<0.2)))


save.image("fig2_table1.Rdata")
