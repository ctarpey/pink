### Standardizing Environmental Data
###   lat, long, temp and precip of populations
### Carolyn Tarpey | October 2015 
### ---------------------------------------

envDat <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/POPINFO.txt", header = TRUE, sep = "\t")
envDat <- as.data.frame(envDat)  
Std_lat<- scale(envDat$latitude)
Std_long<- scale(envDat$Longitude)
envDat$latitude
Std_lat
envDat$Longitude
Std_long


POP_ENV_noheader <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/POP_ENV_noheader.txt", header = FALSE, sep = "\t")
POP_ENV_noheader <- as.data.frame(POP_ENV_noheader) 
colsPOP_ENV_noheader <- c(names(POP_ENV_noheader))
POP_ENV_noheader.dat <- scale(POP_ENV_noheader)
POP_ENV_noheader.dat
write.table(POP_ENV_noheader.dat, file = "POP_ENV_noheader.txt", sep ="\t", row.names = TRUE, col.names = TRUE)


POPs <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/POPs_data.txt", header = FALSE, sep = "\t")
POPs <- as.data.frame(POPs) 
colsPOPs <- c(names(POPs))
POPs.dat <- scale(POPs)
POPs.dat
write.table(POPs.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/POPs.txt", sep ="\t", row.names = TRUE, col.names = TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  O C T O B E R 1 4 

all_env <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/raw_env.txt", header = FALSE, sep = "\t")
all_env <- as.data.frame(all_env) 
colsall_env <- c(names(all_env))
all_env.dat <- scale(all_env)
all_env.dat
write.table(all_env.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/all_env_STAND.txt", sep ="\t", row.names = TRUE, col.names = TRUE)


#########  S T A N D A R D I Z E  T H E  P C A  A X I E S 


PCAs <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/env_PCA_axies.txt", header = FALSE, sep = "\t")
PCAs <- as.data.frame(PCAs) 
colsPCAs <- c(names(PCAs))
PCAs.dat <- scale(PCAs)
PCAs.dat
write.table(PCAs.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/env_PCA_axies_STAND.txt", sep ="\t", row.names = TRUE, col.names = TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  O C T O B E R 3 0

all_env <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/all.txt", header = FALSE, sep = "\t")
all_env <- as.data.frame(all_env) 
colsall_env <- c(names(all_env))
all_env.dat <- scale(all_env)
all_env.dat
write.table(all_env.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/all_stand.txt", sep ="\t", row.names = TRUE, col.names = TRUE)


#########  S T A N D A R D I Z E  T H E  ASIA AND NA seperately

AS_env <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/asial.txt", header = FALSE, sep = "\t")
AS_env <- as.data.frame(AS_env) 
colsAS_env <- c(names(AS_env))
AS_env.dat <- scale(AS_env)
AS_env.dat
write.table(AS_env.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/asial_stand.txt", sep ="\t", row.names = TRUE, col.names = TRUE)

NA_env <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/NA.txt", header = FALSE, sep = "\t")
NA_env <- as.data.frame(NA_env) 
colsNA_env <- c(names(NA_env))
NA_env.dat <- scale(NA_env)
NA_env.dat
write.table(NA_env.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/RestandBayescenv/NA_stand.txt", sep ="\t", row.names = TRUE, col.names = TRUE)





###INDIVIDUAL variables


precip <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/precip.txt", header = FALSE, sep = "\t")
precip <- as.data.frame(precip) 
colsprecip <- c(names(precip))
precip_scaled.dat <- scale(precip)
precip_scaled.dat

tmax <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/tmax.txt", header = FALSE, sep = "\t")
tmax <- as.data.frame(tmax) 
colstmax <- c(names(tmax))
tmax_scaled.dat <- scale(tmax)
tmax_scaled.dat

tmin <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/tmin.txt", header = FALSE, sep = "\t")
tmin <- as.data.frame(tmin) 
colstmin <- c(names(tmin))
tmin_scaled.dat <- scale(tmin)
tmin_scaled.dat

tmean <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/tmean.txt", header = FALSE, sep = "\t")
tmean <- as.data.frame(tmean) 
colstmean <- c(names(tmean))
tmean_scaled.dat <- scale(tmean)
tmean_scaled.dat

bioclim <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/bioclim.txt", header = FALSE, sep = "\t")
bioclim <- as.data.frame(bioclim) 
colsbioclim<- c(names(bioclim))
bioclim_scaled.dat <- scale(bioclim)
bioclim_scaled.dat

write.table(bioclim_scaled.dat, file = "G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/bioclim_stand.txt", sep ="\t", row.names = TRUE, col.names = TRUE)




# check that we get mean of 0 and sd of 1
colMeans(Std_envDat)  # faster version of apply(scaled.dat, 2, mean)
apply(Std_envDat, 2, sd)
