


install.packages("rgdal")
install.packages("raster")

library(rgdal)
library(raster)

setwd("C:\\Users\\Carolyn\\Documents\\GitHub\\pink\\worldclim\\ENVdata")

w = getData(worldclim, download = TRUE, path = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/tmin_30s_bil.zip", var = "tmin") 

w = getData(worldclim, download = TRUE, path = "http://www.worldclim.org", var = "tmin") 

#raster layer
r = raster("tmin12_111.tif.ovr") 
head(r)
cellStats(st,r)

cells <- cellFromRowCol(r, rownr= range(1:20), colnr = range(1:10))

cells <- cellFromXY(r, 57.7, 162.48)

xy =xyFromCell(r, cells)

extract(r, xy)
 ##xy needs to be a matrix of x and y coordinates

cellFromXY(r, c(0,0))