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
setwd("Enter the directory which stores the .RData files from the scripts casex_dony.R")
  
#----------------------------------------------------------

#*************** Case 1 ***************************************
c1<-list(4)
c2<-list(4)
med<-list(4)
truh<-list(4)
donor<-list(4)
for(j in 1:4){
  a=paste('case1_don',j,'_truh.RData',sep='')
  load(a)
  c1[[j]]<-out$c1
  c2[[j]]<-out$c2
  truh[[j]]<-out$teststat
  med[[j]]<-median(out$null.dist)
  donor[[j]]<-j
}
plotdata1<- as.data.frame(cbind(do.call(rbind,c1), do.call(rbind,c2),
                                do.call(rbind,med),do.call(rbind,truh),unlist(donor)))
names(plotdata1)<-c('c1','c2','med','truh','Donor')
plotdata2<-plotdata1[,c(3,5)]
g3<-ggplot()+
  geom_point(data=plotdata1,aes(x=Donor,y=truh,size=2), colour="blue",alpha=1)+
  geom_point(data=plotdata2,aes(x=Donor,y=med,size=2),shape=8,color="red")+
  geom_errorbar(data=plotdata1,width=.1, aes(x=Donor,ymin=c1, ymax=c2),size=1, colour="red")+ylim(1,8)+
  xlab('Donor')+ylab("Null Distribution")+theme_bw()+
  ggtitle("CASE 1: Uninfected versus Nef rich HIV Infected (d=35)")+
  theme(legend.position="none",
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5))

#*************** Case 2 ***************************************
c1<-list(4)
c2<-list(4)
med<-list(4)
truh<-list(4)
donor<-list(4)
for(j in 1:4){
  a=paste('case2_don',j,'_truh_new.RData',sep='')
  load(a)
  c1[[j]]<-out$c1
  c2[[j]]<-out$c2
  truh[[j]]<-out$teststat
  med[[j]]<-median(out$null.dist)
  donor[[j]]<-j
}
plotdata1<- as.data.frame(cbind(do.call(rbind,c1), do.call(rbind,c2),
                                do.call(rbind,med),do.call(rbind,truh),unlist(donor)))
names(plotdata1)<-c('c1','c2','med','truh','Donor')
plotdata2<-plotdata1[,c(3,5)]
g5<-ggplot()+
  geom_point(data=plotdata1,aes(x=Donor,y=truh,size=2), colour="blue",alpha=1)+
  geom_point(data=plotdata2,aes(x=Donor,y=med,size=2),shape=8,color="red")+
  geom_errorbar(data=plotdata1,width=.1, aes(x=Donor,ymin=c1, ymax=c2),size=1, colour="red")+ylim(1,8)+
  xlab('Donor')+ylab("Null Distribution")+theme_bw()+
  ggtitle("CASE 2: Uninfected versus Nef rich HIV Infected (d=31)")+
  theme(legend.position="none",
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5))

grid.arrange(g3,g5,ncol=2)

#*************** Case 3 ***************************************
c1<-list(4)
c2<-list(4)
med<-list(4)
truh<-list(4)
donor<-list(4)
for(j in 1:4){
  a=paste('case3_don',j,'_truh.RData',sep='')
  load(a)
  c1[[j]]<-out$c1
  c2[[j]]<-out$c2
  truh[[j]]<-out$teststat
  med[[j]]<-median(out$null.dist)
  donor[[j]]<-j
}
plotdata1<- as.data.frame(cbind(do.call(rbind,c1), do.call(rbind,c2),
                                do.call(rbind,med),do.call(rbind,truh),unlist(donor)))
names(plotdata1)<-c('c1','c2','med','truh','Donor')
plotdata2<-plotdata1[,c(3,5)]
g4<-ggplot()+
  geom_point(data=plotdata1,aes(x=Donor,y=truh,size=2), colour="blue",alpha=1)+
  geom_point(data=plotdata2,aes(x=Donor,y=med,size=2),shape=8,color="red")+
  geom_errorbar(data=plotdata1,width=.1, aes(x=Donor,ymin=c1, ymax=c2),size=1, colour="red")+ylim(1,8)+
  xlab('Donor')+ylab("Null Distribution")+theme_bw()+
  ggtitle("CASE 3: Uninfected versus Nef deficient HIV Infected (d=35)")+
  theme(legend.position="none",
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5))

#*************** Case 4 ***************************************
c1<-list(4)
c2<-list(4)
med<-list(4)
truh<-list(4)
donor<-list(4)
for(j in 1:4){
  a=paste('case4_don',j,'_truh.RData',sep='')
  load(a)
  c1[[j]]<-out$c1
  c2[[j]]<-out$c2
  truh[[j]]<-out$teststat
  med[[j]]<-median(out$null.dist)
  donor[[j]]<-j
}
plotdata1<- as.data.frame(cbind(do.call(rbind,c1), do.call(rbind,c2),
                                do.call(rbind,med),do.call(rbind,truh),unlist(donor)))
names(plotdata1)<-c('c1','c2','med','truh','Donor')
plotdata2<-plotdata1[,c(3,5)]
g6<-ggplot()+
  geom_point(data=plotdata1,aes(x=Donor,y=truh,size=2), colour="blue",alpha=1)+
  geom_point(data=plotdata2,aes(x=Donor,y=med,size=2),shape=8,color="red")+
  geom_errorbar(data=plotdata1,width=.1, aes(x=Donor,ymin=c1, ymax=c2),size=1, colour="red")+ylim(1,8)+
  xlab('Donor')+ylab("Null Distribution")+theme_bw()+
  ggtitle("CASE 4: Uninfected versus Nef deficient HIV Infected (d=31)")+
  theme(legend.position="none",
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.title = element_text(hjust = 0.5))

grid.arrange(g4,g6,ncol=2)



