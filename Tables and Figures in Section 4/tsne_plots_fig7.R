
library(Rtsne)
library(ggplot2)
library(gridExtra)

source('lib.R')
datalocation<- "Infection Data" # Please enter the data directory here

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
tsne_plot<-function(datalocation,datafilename1,datafilename2,drop_list,donor,case){
  
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
  UN<-D_UN[,index];
  V<- D_V[,index];
  
  UN<-D_UN[,2:39];
  V<- D_V[,2:39];
  index1<- which(names(UN) %in% drop_list)
  index2<- which(names(V) %in% drop_list)
  V.sub<-data.matrix(V[,-index2]);
  UN.sub<-data.matrix(UN[,-index1]);
  n.u<-dim(UN.sub)[1]
  n.v<-dim(V.sub)[1]
  
  UV<-as.matrix(rbind(UN.sub,V.sub))
  
  set.seed(42)
  tsne_out<-Rtsne(UV)
  dat1<-as.data.frame(tsne_out$Y)
  dat1$Case<-as.factor(c(rep('Uninfected',n.u),rep('Infected',n.v)))
  g<-ggplot() +
    geom_point(data=dat1,aes(x=V1,y=V2,color=Case,shape=Case),size=0.7)+
    xlab(paste("Case",case,", Donor ", donor))+ylab("")+theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  return(g)
  
}

#-------------------------------------------------------
#1. Case 1 donors 1-4
#----------------------------------------------------------
datafilename1<- "1.D1UN_Uninfected CD4s.txt"
datafilename2<- "2.D1F4_HSA(hi).txt"
drop_list<- c('CD3','CD19','HSA')
g3_1<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,1,1)

datafilename1<- "4.D2UN_Uninfected CD4s.txt"
datafilename2<- "5.D2F4_HSA(hi).txt"
g3_2<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,2,1)

datafilename1<- "7.D3UN_Uninfected CD4s.txt"
datafilename2<- "8.D3F4_HSA(hi).txt"
g3_3<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,3,1)

datafilename1<- "10.D4UN_Uninfected CD4s.txt"
datafilename2<- "11.D4F4_HSA(hi).txt"
g3_4<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,4,1)

#--------------------------------------------------------------
#2. Case 2 donors 1-4
#----------------------------------------------------------
datafilename1<- "1.D1UN_Uninfected CD4s.txt"
datafilename2<- "2.D1F4_HSA(hi).txt"
drop_list<- c('CD3','CD19','HSA','CD4','CCR5','CD28','CD62L')
g5_1<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,1,2)

datafilename1<- "4.D2UN_Uninfected CD4s.txt"
datafilename2<- "5.D2F4_HSA(hi).txt"
g5_2<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,2,2)

datafilename1<- "7.D3UN_Uninfected CD4s.txt"
datafilename2<- "8.D3F4_HSA(hi).txt"
g5_3<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,3,2)

datafilename1<- "10.D4UN_Uninfected CD4s.txt"
datafilename2<- "11.D4F4_HSA(hi).txt"
g5_4<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,4,2)
#--------------------------------------------------------------

#3. Case 3 and donors 1-4
#----------------------------------------------------------
datafilename1<- "1.D1UN_Uninfected CD4s.txt"
datafilename2<- "3.D1dNef_HSA(hi).txt"
drop_list<- c('CD3','CD19','HSA')
g4_1<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,1,3)

datafilename1<- "4.D2UN_Uninfected CD4s.txt"
datafilename2<- "6.D2dNef_HSA(hi).txt"
g4_2<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,2,3)

datafilename1<- "7.D3UN_Uninfected CD4s.txt"
datafilename2<- "9.D3dNef_HSA(hi).txt"
g4_3<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,3,3)

datafilename1<- "10.D4UN_Uninfected CD4s.txt"
datafilename2<- "12.D4dNef_HSA(hi).txt"
g4_4<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,4,3)
#--------------------------------------------------------------
#4. Case 4 and donors 1-4
#----------------------------------------------------------
datafilename1<- "1.D1UN_Uninfected CD4s.txt"
datafilename2<- "3.D1dNef_HSA(hi).txt"
drop_list<- c('CD3','CD19','HSA','CD4','CCR5','CD28','CD62L')
g6_1<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,1,4)

datafilename1<- "4.D2UN_Uninfected CD4s.txt"
datafilename2<- "6.D2dNef_HSA(hi).txt"
g6_2<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,2,4)

datafilename1<- "7.D3UN_Uninfected CD4s.txt"
datafilename2<- "9.D3dNef_HSA(hi).txt"
g6_3<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,3,4)

datafilename1<- "10.D4UN_Uninfected CD4s.txt"
datafilename2<- "12.D4dNef_HSA(hi).txt"
g6_4<-tsne_plot(datalocation,datafilename1,datafilename2,drop_list,4,4)
#--------------------------------------------------------------

#### build the plot #########

ggarrange(g3_1, g3_2, g3_3, g3_4, g5_1, g5_2, g5_3, g5_4,
          g4_1, g4_2, g4_3, g4_4, g6_1, g6_2, g6_3, g6_4,
          ncol=4, nrow=4, common.legend = TRUE, legend="top")

save.image("tsne_plot.RData")

