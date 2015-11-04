### PCA of Environmental Data
###   choosing uncorrelated environmental variables to continue analysis
### Carolyn Tarpey | October 2015 
### ---------------------------------------


install.packages("FactoMineR")
library(FactoMineR)
library(ggplot2)

popnames <-  read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/POPNAMES.txt", header = FALSE, sep = "\t")

env_data <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_HORZ.txt", header = TRUE, sep = "\t")
env_data  <- as.data.frame(env_data)  

names(env_data)
#View(env_data)

pca_env <- prcomp(env_data)

pca_env$sdev

#matrix whose column contain the eigenvectors
head(pca_env$rotation)

head(pca_env$x)

summary(pca_env)
biplot(pca_env)

# produce a scree plot 
plot(pca_env)

#create data frame with scores
scores <- as.data.frame(pca_env$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables") + theme_classic()



circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data, pca_env$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), 
                    y1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Environmental Correlations") + theme_classic()

write.table(pca_env$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores.txt", sep ="\t", row.names = TRUE, col.names = TRUE )








################################################### Mix no bio

env_data_mix_noBio <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_mix_noBio.txt", header = TRUE, sep = "\t")
env_data_mix_noBio  <- as.data.frame(env_data_mix_noBio)  

names(env_data_mix_noBio)
#View(env_data)

pca_env_mix_noBio <- prcomp(env_data_mix_noBio)

pca_env_mix_noBio$sdev

#matrix whose column contain the eigenvectors
head(pca_env_mix_noBio$rotation)

head(pca_env_mix_noBio$x)

summary(pca_env_mix_noBio)
biplot(pca_env_mix_noBio)

# produce a scree plot 
plot(pca_env_mix_noBio)

#create data frame with scores
scores_mix_noBio <- as.data.frame(pca_env_mix_noBio$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_mix_noBio, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables") + theme_classic()



circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_mix_noBio, pca_env_mix_noBio$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Environmental Correlations") + theme_classic()

write.table(pca_env_mix_noBio$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_mix_noBio.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_mix_noBio$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_mix_noBio.txt", sep ="\t", row.names = TRUE, col.names = TRUE )




 


################################################### Mix 

env_data_mix <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_mix.txt", header = TRUE, sep = "\t")
env_data_mix  <- as.data.frame(env_data_mix)  

names(env_data_mix)
#View(env_data)

pca_env_mix <- prcomp(env_data_mix)

pca_env_mix$sdev

#matrix whose column contain the eigenvectors
head(pca_env_mix$rotation)

head(pca_env_mix$x)

summary(pca_env_mix)
biplot(pca_env_mix)

# produce a scree plot 
plot(pca_env_mix)

#create data frame with scores
scores_mix <- as.data.frame(pca_env_mix$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_mix, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental Variables") + theme_classic()



circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_mix, pca_env_mix$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Environmental Correlations") + theme_classic()

write.table(pca_env_mix$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_mix.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_mix$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_mix.txt", sep ="\t", row.names = TRUE, col.names = TRUE )


















########################  Temp

env_data_temp <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_Temp.txt", header = TRUE, sep = "\t")
env_data_temp  <- as.data.frame(env_data_temp )  

names(env_data_temp)
#View(env_data)

pca_env_temp <- prcomp(env_data_temp)

pca_env_temp$sdev
head(pca_env_temp$rotation)
head(pca_env_temp$x)

summary(pca_env_temp)
biplot(pca_env_temp)

# produce a scree plot 
plot(pca_env_temp)

#create data frame with scores
scores_temp <- as.data.frame(pca_env_temp$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_temp, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Temperature") + theme_classic()

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_temp, pca_env_temp$x))

head(correlations)


# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Correlations: Temperature") + theme_classic()


write.table(pca_env_temp$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_temp.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_temp$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_temp.txt", sep ="\t", row.names = TRUE, col.names = TRUE )


########################  Precip

env_data_precip <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_Precip.txt", header = TRUE, sep = "\t")
env_data_precip  <- as.data.frame(env_data_precip )  

names(env_data_precip)
#View(env_data)

pca_env_precip <- prcomp(env_data_precip)

pca_env_precip$sdev
head(pca_env_precip$rotation)
head(pca_env_precip$x)

summary(pca_env_precip)
biplot(pca_env_precip)

#create data frame with scores
scores_precip <- as.data.frame(pca_env_precip$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_precip, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Precipitation") + theme_classic()

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_precip, pca_env_precip$x))

head(correlations)


# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Correlations: Precipitation") + theme_classic()


write.table(pca_env_precip$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_precip.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_precip$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_precip.txt", sep ="\t", row.names = TRUE, col.names = TRUE )



########################  Bioclim


env_data_Bioclim <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_Bioclim.txt", header = TRUE, sep = "\t")
env_data_Bioclim  <- as.data.frame(env_data_Bioclim )  

names(env_data_Bioclim)
#View(env_data)

pca_env_Bioclim <- prcomp(env_data_Bioclim) 

pca_env_Bioclim$sdev
head(pca_env_Bioclim$rotation)
head(pca_env_Bioclim$x)

summary(pca_env_Bioclim)
biplot(pca_env_Bioclim)

#create data frame with scores
scores_Bioclim <- as.data.frame(pca_env_Bioclim$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_Bioclim, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Bioclim") + theme_classic()

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_Bioclim, pca_env_Bioclim$x))

head(correlations)


# data frame with arrows coordinates
arrows = data.frame(x1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of correlations") + theme_classic()


write.table(pca_env_Bioclim$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_Bioclim.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_Bioclim$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_Bioclim.txt", sep ="\t", row.names = TRUE, col.names = TRUE )







######################## NO Bioclim


env_data_noBio <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_noBio.txt", header = TRUE, sep = "\t")
env_data_noBio  <- as.data.frame(env_data_noBio )  

names(env_data_noBio)
#View(env_data)

pca_env_noBio <- prcomp(env_data_noBio) 

pca_env_noBio$sdev
head(pca_env_noBio$rotation)
head(pca_env_noBio$x)

summary(pca_env_noBio)
biplot(pca_env_noBio)

#create data frame with scores
scores_noBio <- as.data.frame(pca_env_noBio$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_noBio, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Environmental variables without Bioclim") + theme_classic()

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_noBio, pca_env_noBio$x))

head(correlations)


# data frame with arrows coordinates
arrows = data.frame(x1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 
                    y1 = c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Correlations: Not including Bioclim") + theme_classic()


write.table(pca_env_noBio$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_noBio.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_noBio$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_noBio.txt", sep ="\t", row.names = TRUE, col.names = TRUE )


########################  Tmean

env_data_Tmean <- read.table ("G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_4_pca_Tmean.txt", header = TRUE, sep = "\t")
env_data_Tmean  <- as.data.frame(env_data_Tmean )  

names(env_data_Tmean)
#View(env_data)

pca_env_Tmean <- prcomp(env_data_Tmean)

pca_env_Tmean$sdev
head(pca_env_Tmean$rotation)
head(pca_env_Tmean$x)

summary(pca_env_Tmean)
biplot(pca_env_Tmean)

#create data frame with scores
scores_Tmean <- as.data.frame(pca_env_Tmean$x)

#plot the observations ###THIS IS THE NAMES OF THE POPS
ggplot(data = scores_Tmean, aes(x= PC1, y = PC2, label = popnames[,1])) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 10) +
  ggtitle("PCA plot of Mean Temperature") + theme_classic()

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(env_data_Tmean, pca_env_Tmean$x))

head(correlations)


# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ), 
                    y1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ),
                    x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+ 
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs",  y = "pc2 axis") + 
  ggtitle("Circle of Correlations: Mean Temperature") + theme_classic()


write.table(pca_env_Tmean$rotation, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_rotation_Tmean.txt", sep ="\t", row.names = TRUE, col.names = TRUE )
write.table(pca_env_Tmean$x, file ="G:/Analysis/Pop_analysis/Populations_b3_may/Env_data_pca_scores_Tmean.txt", sep ="\t", row.names = TRUE, col.names = TRUE )











###USING ANOTHER METHOD
pcaFM <- PCA(env_data, graph = FALSE)
#matrix with eigenvalues
pcaFM$eig

#create data framw with scores
scoresFM <- as.data.frame(pcaFM$var)
names(pcaFM$var)

#plot the observations
ggplot(data = scoresFM, aes(x= PC1, y = PC2, label = rownames(scoresFM))) +
  geom_hline(yintercept = 0, color = "gray65") +
  geom_vline(xintercept = 0, color = "gray65") +
  geom_text(color = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Environmental Variables")

