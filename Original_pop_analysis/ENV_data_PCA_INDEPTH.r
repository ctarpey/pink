### PCA of Environmental Data
###   More in depth: lat long Temp and precip
### Carolyn Tarpey | October 2015 
### ---------------------------------------


library(FactoMineR)
library(ggplot2)

popnames <-  read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/POPNAMES.txt", header = FALSE, sep = "\t")

Precip <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Precip.txt", header = TRUE, sep = "\t")
Precip_lat <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Precip_lat.txt", header = TRUE, sep = "\t")
Temp <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Temp.txt", header = TRUE, sep = "\t")
Temp_lat <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Temp_lat.txt", header = TRUE, sep = "\t")
Precip_temp <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Precip_temp.txt", header = TRUE, sep = "\t")
Precip_temp_lat <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/R_output/Precip_temp_lat.txt", header = TRUE, sep = "\t")
Bio_sm <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/Bio_sm.txt", header = TRUE, sep = "\t")
Bio_lat_sm <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/Bio_lat_sm.txt", header = TRUE, sep = "\t")


Precip <- as.data.frame(Precip) 
Precip_lat <- as.data.frame(Precip_lat) 
Temp <- as.data.frame(Temp) 
Temp_lat <- as.data.frame(Temp_lat) 
Precip_temp <- as.data.frame(Precip_temp) 
Precip_temp_lat <- as.data.frame(Precip_temp_lat) 
Bio_sm <- as.data.frame(Bio_sm) 
Bio_lat_sm <- as.data.frame(Bio_lat_sm) 


names(Precip) 
names(Precip_lat) 
names(Temp) 
names(Temp_lat) 
names(Precip_temp) 
names(Precip_temp_lat) 
names(Bio_sm) 
names(Bio_lat_sm) 

pca_Precip <- prcomp(Precip) 
pca_Precip_lat <- prcomp(Precip_lat)
pca_Temp <- prcomp(Temp) 
pca_Temp_lat <- prcomp(Temp_lat)  
pca_Precip_temp <- prcomp(Precip_temp) 
pca_Precip_temp_lat <- prcomp(Precip_temp_lat)
pca_Bio_sm <- prcomp(Bio_sm) 
pca_Bio_lat_sm <- prcomp(Bio_lat_sm)


pca_Precip$sdev 
pca_Precip_lat$sdev
pca_Temp$sdev
pca_Temp_lat$sdev
pca_Precip_temp$sdev
pca_Precip_temp_lat$sdev
pca_Bio_sm$sdev
pca_Bio_lat_sm$sdev

rotation_p <- pca_Precip$rotation 
rotation_p_l <- pca_Precip_lat$rotation
rotation_t <- pca_Temp$rotation
rotation_t_l <- pca_Temp_lat$rotation
rotation_p_t <- pca_Precip_temp$rotation
rotation_p_t_l <- pca_Precip_temp_lat$rotation
rotation_b <- pca_Bio_sm$rotation
rotation_b_l <- pca_Bio_lat_sm$rotation

#matrix whose column contain the eigenvectors
#head(pca_env$rotation)
#head(pca_env$x)

# summary stats of PCA
summary(pca_Precip)
summary(pca_Precip_lat)
summary(pca_Temp)
summary(pca_Temp_lat)
summary(pca_Precip_temp)
summary(pca_Precip_temp_lat)
summary(pca_Bio_sm)
summary(pca_Bio_lat_sm)

# produce a biplot plot
biplot(pca_Precip)
biplot(pca_Precip_lat)
biplot(pca_Temp)
biplot(pca_Temp_lat)
biplot(pca_Precip_temp)
biplot(pca_Precip_temp_lat)
biplot(pca_Bio_sm)
biplot(pca_Bio_lat_sm)

# produce a scree plot
plot(pca_Precip)
plot(pca_Precip_lat)
plot(pca_Temp)
plot(pca_Temp_lat)
plot(pca_Precip_temp)
plot(pca_Precip_temp_lat)
plot(pca_Bio_sm)
plot(pca_Bio_lat_sm)

#create data frame with scores
scores_p <- as.data.frame(pca_Precip$x)
scores_p_l <- as.data.frame(pca_Precip_lat$x)
scores_t <- as.data.frame(pca_Temp$x)
scores_t_l <- as.data.frame(pca_Temp_lat$x)
scores_p_t <- as.data.frame(pca_Precip_temp$x)
scores_p_t_l <- as.data.frame(pca_Precip_temp_lat$x)
scores_b <- as.data.frame(pca_Bio_sm$x)
scores_b_l <- as.data.frame(pca_Bio_lat_sm$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_p, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Precip") + theme_classic()

ggplot(data = scores_p_l, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Precip Lat long") + theme_classic()

ggplot(data = scores_t, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Temp") + theme_classic()

ggplot(data = scores_t_l, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Temp lat long") + theme_classic()

ggplot(data = scores_p_t, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Precip Temp ") + theme_classic()

ggplot(data = scores_p_t_l, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Precip Temp Lat Long ") + theme_classic()

ggplot(data = scores_b, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Bio ") + theme_classic()

ggplot(data = scores_b_l, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables Bio Lat Long ") + theme_classic()



#plot the evnironmental variables

#make the loadings into a dataframe
Dfrotation_p <- as.data.frame(rotation_p)
Dfrotation_p_l <- as.data.frame(rotation_p_l)
DFrotation_t <- as.data.frame(rotation_t)
DFrotation_t_l <- as.data.frame(rotation_t_l)
DFrotation_p_t <- as.data.frame(rotation_p_t)
DFrotation_p_t_l <- as.data.frame(rotation_p_t_l)
DFrotation_b <- as.data.frame(rotation_b)
DFrotation_b_l <- as.data.frame(rotation_b_l)


ggplot(data = Dfrotation_p, aes(x= PC1, y = PC2, label = rownames(Dfrotation_p) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip ") + theme_classic()

ggplot(data = Dfrotation_p_l, aes(x= PC1, y = PC2, label = rownames(Dfrotation_p_l) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip ") + theme_classic()

ggplot(data = DFrotation_t, aes(x= PC1, y = PC2, label = rownames(DFrotation_t) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip ") + theme_classic()

ggplot(data = DFrotation_t_l, aes(x= PC1, y = PC2, label = rownames(DFrotation_t_l) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip ") + theme_classic()

ggplot(data = DFrotation_p_t, aes(x= PC1, y = PC2, label = rownames(DFrotation_p_t) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip ") + theme_classic()

ggplot(data = DFrotation_p_t_l, aes(x= PC1, y = PC2, label = rownames(DFrotation_p_t_l) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Precip Temp lat long ") + theme_classic()

ggplot(data = DFrotation_b, aes(x= PC1, y = PC2, label = rownames(DFrotation_b) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Bio ") + theme_classic()

ggplot(data = DFrotation_b_l, aes(x= PC1, y = PC2, label = rownames(DFrotation_b_l) )) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables Bio lat long ") + theme_classic()


#loading plot
plot(pca_Precip$rotation[,1], pca_Precip$rotation[,2],
     xlim=c(-1,1), ylim=c(-1,1),
     main='Loadings for PC1 vs. PC2 Precipitation')
text(rotation_p[,1]*1.55,rotation_p[,2], 
     rownames(rotation_p), col="red", cex=0.7)

plot(pca_Temp$rotation[,1], pca_Temp$rotation[,2],
     xlim=c(-1,1), ylim=c(-1,1),
     main='Loadings for PC1 vs. PC2 Temperature')
text(rotation_t[,1]*1.55,rotation_t[,2], 
     rownames(rotation_t), col="red", cex=0.7)

plot(pca_Precip_temp$rotation[,1], pca_Precip_temp$rotation[,2],
     xlim=c(-.8,.5), ylim=c(-.5,.5),
     main='Loadings for PC1 vs. PC2 Precipitation and Temperature')
text(rotation_p_t[,1]*1.8,rotation_p_t[,2], 
     rownames(rotation_p_t), col="red", cex=0.7)

plot(pca_Precip_lat$rotation[,1], pca_Precip_lat$rotation[,2],
     xlim=c(-1,1), ylim=c(-1,1),
     main='Loadings for PC1 vs. PC2 Precipitation and Latitude')
text(rotation_p_l[,1]*1.55,rotation_p_l[,2], 
     rownames(rotation_p_l), col="red", cex=0.7)

plot(pca_Temp_lat$rotation[,1], pca_Temp_lat$rotation[,2],
     xlim=c(-1,1), ylim=c(-1,1),
     main='Loadings for PC1 vs. PC2 Temperature and latitude')
text(rotation_t_l[,1]*1.55,rotation_t_l[,2], 
     rownames(rotation_t_l), col="red", cex=0.7)

plot(pca_Precip_temp_lat$rotation[,1], pca_Precip_temp_lat$rotation[,2],
     xlim=c(-.8,.5), ylim=c(-.5,.5),
     main='Loadings for PC1 vs. PC2 Precipitation Temperature and Latitude')
text(rotation_p_t_l[,1]*1.8,rotation_p_t_l[,2], 
     rownames(rotation_p_t_l), col="red", cex=0.7)


#correlation loading plot

correlationloadings_precip <- cor(Precip, pca_Precip$x)
plot(correlationloadings_precip[,1], correlationloadings_precip[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Precipitation')
text(correlationloadings_precip[,1]*1.3,correlationloadings_precip[,2], 
     rownames(rotation_p), col="red", cex=0.7)

correlationloadings_temp <- cor(Temp, pca_Temp$x)
plot(correlationloadings_temp[,1], correlationloadings_temp[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Temperature')
text(correlationloadings_temp[,1]*1.3,correlationloadings_temp[,2], 
     rownames(rotation_t), col="red", cex=0.7)

correlationloadings_Precip_temp <- cor(Precip_temp, pca_Precip_temp$x)
plot(correlationloadings_Precip_temp[,1], correlationloadings_Precip_temp[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Precipitation and Temperature')
text(correlationloadings_Precip_temp[,1]*1.3,correlationloadings_Precip_temp[,2], 
     rownames(rotation_p_t), col="red", cex=0.7)

correlationloadings_Precip_lat <- cor(Precip_lat, pca_Precip_lat$x)
plot(correlationloadings_Precip_lat[,1], correlationloadings_Precip_lat[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Precipitation Temperature and Latitude')
text(correlationloadings_Precip_lat[,1]*1.3,correlationloadings_Precip_lat[,2], 
     rownames(rotation_p_l), col="red", cex=0.7)

correlationloadings_Temp_lat <- cor(Temp_lat, pca_Temp_lat$x)
plot(correlationloadings_Temp_lat[,1], correlationloadings_Temp_lat[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Precipitation Temperature and Latitude')
text(correlationloadings_Temp_lat[,1]*1.3,correlationloadings_Temp_lat[,2], 
     rownames(rotation_t_l), col="red", cex=0.7)

correlationloadings_Precip_temp_lat <- cor(Precip_temp_lat, pca_Precip_temp_lat$x)
plot(correlationloadings_Precip_temp_lat[,1], correlationloadings_Precip_temp_lat[,2],
     xlim=c(-1.5,1), ylim=c(-1,1),
     main='Correlation Loadings for PC1 vs. PC2 Precipitation Temperature and Latitude')
text(correlationloadings_Precip_temp_lat[,1]*1.3,correlationloadings_Precip_temp_lat[,2], 
     rownames(rotation_p_t_l), col="red", cex=0.7)


#########################################################
# this opens a new graph window before creating a graph so you dont overwrite the other graph 

windows(height=7, width=7)

plot(scores_p[,1], scores_p[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 3), ylim =c(-4.5, 3))
arrows(0,0, rotation_p[,1]*8,rotation_p[,2]*8, length=0.1, angle=20, col="red")
text(rotation_p[,1]*8*1.3,rotation_p[,2]*8*1.02, 
     rownames(rotation_p), col="red", cex=0.7)


plot(scores_p_l[,1], scores_p_l[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 3), ylim =c(-4.5, 5))
arrows(0,0, rotation_p_l[,1]*8,rotation_p_l[,2]*8, length=0.1, angle=20, col="red")
text(rotation_p_l[,1]*8*1.3,rotation_p_l[,2]*8*1.02, 
     rownames(rotation_p_l), col="red", cex=0.7)

plot(scores_t[,1], scores_t[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 3), ylim =c(-4.5, 4))
arrows(0,0, rotation_t[,1]*8,rotation_t[,2]*8, length=0.1, angle=20, col="red")
text(rotation_t[,1]*8*1.3,rotation_t[,2]*8*1.02, 
     rownames(rotation_t), col="red", cex=0.7)

plot(scores_t_l[,1], scores_t_l[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 4), ylim =c(-4.5, 5))
arrows(0,0, rotation_t_l[,1]*8,rotation_t_l[,2]*8, length=0.1, angle=20, col="red")
text(rotation_t_l[,1]*8*1.3,rotation_t_l[,2]*8*1.02, 
     rownames(rotation_t_l), col="red", cex=0.7)

plot(scores_p_t[,1], scores_p_t[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 2.5), ylim =c(-2.5, 3))
arrows(0,0, rotation_p_t[,1]*8,rotation_p_t[,2]*8, length=0.1, angle=20, col="red")
text(rotation_p_t[,1]*8*1.3,rotation_p_t[,2]*8*1.02, 
     rownames(rotation_p_t), col="red", cex=0.7)

plot(scores_p_t_l[,1], scores_p_t_l[,2], xlab="PCA 1", ylab="PCA 2", 
     type="n", xlim=c(-4.5, 4), ylim =c(-2.5, 2.5))
arrows(0,0, rotation_p_t_l[,1]*8,rotation_p_t_l[,2]*8, length=0.1, angle=20, col="red")
text(rotation_p_t_l[,1]*8*1.3,rotation_p_t_l[,2]*8*1.02, 
     rownames(rotation_p_t_l), col="red", cex=0.7)


#########################################################


write.table(pca_Precip$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Precip$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Precip_lat$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_lat_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Precip_lat$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_lat_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Temp$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Temp_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Temp$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Temp_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Temp_lat$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Temp_lat_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Temp_lat$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Temp_lat_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Precip_temp_lat$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_temp_lat_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Precip_temp_lat$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_temp_lat_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Precip_temp$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_temp_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Precip_temp$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/pca_Precip_temp_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Bio_sm$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/pca_Precip_temp_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Bio_sm$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/pca_Precip_temp_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )

write.table(pca_Bio_lat_sm$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/pca_Bio_lat_sm_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_Bio_lat_sm$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/EnvData/VariableTxtFiles/pca_Bio_lat_sm_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )


