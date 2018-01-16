### Comparing LepMap and MSTMap outputs
###   MSTMap run by family, LepMap of all families together
### Carolyn Tarpey | August 2015 
### ---------------------------------------

library(lattice)

data<-read.delim("G:/Analysis/Mapping/AllHaps/ComparingLepMapMSTmapLG13.txt",header=TRUE)
data_paralogs <-read.delim("G:/Analysis/Mapping/AllHaps/ComparingLepMapMSTmapLG13paralogs.txt", header=TRUE)

data_paralogs_fam01 <-read.delim("G:/Analysis/Mapping/AllHaps/ComparingLepMapMSTmapLG13_fam01.txt", header=TRUE)

names(data)

xyplot(data$position~data$X05H_5Pos,group=data$X05H_5)
xyplot(data$position~data$X01H_4Pos,group=data$X01H_4)
xyplot(data$position~data$X108_80GRPos,group=data$X108_80GR)
xyplot(data$position~data$X110_65GRPos,group=data$X110_65GR)


names(data_paralogs)

xyplot(data_paralogs$position~data_paralogs$X05H_5Pos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X01H_4Pos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X108_80GRPos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X110_65GRPos,group=data_paralogs$paralog)


names(data_paralogs_fam01)

xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X05H_5Pos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X01H_4Pos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X108_80GRPos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X110_65GRPos,group=data_paralogs_fam01$Just_01)






data_B<-read.delim("G:/Analysis/Mapping/AllHaps/LepMap/ComparingLepNo01to01MSTMap.txt",header=TRUE)

names(data_B)

xyplot(data_B$position~data$X05H_5Pos,group=data$X05H_5)



names(data_paralogs)

xyplot(data_paralogs$position~data_paralogs$X05H_5Pos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X01H_4Pos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X108_80GRPos,group=data_paralogs$paralog)
xyplot(data_paralogs$position~data_paralogs$X110_65GRPos,group=data_paralogs$paralog)


names(data_paralogs_fam01)

xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X05H_5Pos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X01H_4Pos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X108_80GRPos,group=data_paralogs_fam01$Just_01)
xyplot(data_paralogs_fam01$position~data_paralogs_fam01$X110_65GRPos,group=data_paralogs_fam01$Just_01)
