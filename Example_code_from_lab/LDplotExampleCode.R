library(ggplot2)
LDdata<-read.delim("E:/CHUM/WA_chumCompare/batch_4_genepop_omy28loci_plink.ld",sep="")
ggplot()+geom_point(data=LDdata,aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()
ggplot()+geom_point(data=LDdata[LDdata$R2>=0.3,],aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()
