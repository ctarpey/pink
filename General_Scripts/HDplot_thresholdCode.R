thresh_H<-0.65
thresh_H_divDup<-0.88
thresh_D<-5
thresh_Dneg<--5

#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  H<-as.numeric(data[7])
  D<-as.numeric(data[9])
  if(is.na(D)){
    paralog<-NA
  }else if(H>=thresh_H_divDup){
    paralog<-2
  }else if(H>=thresh_H|D>thresh_D|D<thresh_Dneg){
    paralog<-1
  }else{
    paralog<-0
  }
  return(paralog)
}
TogNec_HDplotData$paralog<-apply(TogNec_HDplotData,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg,thresh_H_divDup)


ggplot()+geom_point(data=FigureS1data,aes(x=het_perc,y=z,color=Locus),alpha=0.6)+theme_bw()+
  xlab("H")+ylab("D")+scale_color_manual(values=c("blue","red","dark green"))+
  theme(axis.text=element_text(size=11),axis.title=element_text(size=14,face="bold"))+
  labs(color="Locus Type")
dev.off()