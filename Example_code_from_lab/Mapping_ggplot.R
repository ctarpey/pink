### Creating Geographical Maps in GGPlot
###     
###    This is Eleni Petrou's code for plotting Geographical maps in R
###    
### Eleni Petrou and Carolyn Tarpey | June 2018
## ---------------------------------------

install.packages("tidyverse")
devtools::install_github("tidyverse/lubridate")
install.packages("lubridate")
remove.packages("lubridate")
Features
# Libraries
library(tidyverse)
library(maps)
library(maptools)
library(ggrepel)
library(viridis)
library(ggplot2)
library(pals)
library(grid)
sessionInfo()
.libPaths()
# Get the world polygon and extract USA and Canada

USA <- map_data("world") %>% 
  filter(region=="USA")

Canada <- map_data("world") %>% 
  filter(region=="Canada")


# Read in your data frame with longitude, latitude, and other metadata.
my_data <- read.delim("MAPPING_Herring_PopulationStructure_SamplingLocations.txt")
head(my_data)

# subset the data to get WA samples only
WA_data <- my_data[my_data$Region=="Washington",]

# subset the data to get BC samples only
BC_data <- my_data[my_data$Region=="British Columbia",]



########################################################################################################
# THESE ARE THE FUNCTIONS YOU NEED TO MAKE THE SCALE BAR
#The code in the chunk below was written by Ewen Gallic (bless his heart for being so sharing).

# Result #
#--------#
# Return a list whose elements are :
# 	- rectangle : a data.frame containing the coordinates to draw the first rectangle ;
# 	- rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
# 	- legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_scale_bar <- function(lon,lat,distance_lon,distance_lat,distance_legend, dist_units = "km"){
  # First rectangle
  bottom_right <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon, dist.units = dist_units, model = "WGS84")
  
  topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_lat, dist.units = dist_units, model = "WGS84")
  rectangle <- cbind(lon=c(lon, lon, bottom_right[1,"long"], bottom_right[1,"long"], lon),
                     lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
  rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
  
  # Second rectangle t right of the first rectangle
  bottom_right2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distance_lon*2, dist.units = dist_units, model = "WGS84")
  rectangle2 <- cbind(lon = c(bottom_right[1,"long"], bottom_right[1,"long"], bottom_right2[1,"long"], bottom_right2[1,"long"], bottom_right[1,"long"]),
                      lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
  rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
  
  # Now let's deal with the text
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance_legend, dist.units = dist_units, model = "WGS84")
  on_top2 <- on_top3 <- on_top
  on_top2[1,"long"] <- bottom_right[1,"long"]
  on_top3[1,"long"] <- bottom_right2[1,"long"]
  
  legend <- rbind(on_top, on_top2, on_top3)
  legend <- data.frame(cbind(legend, text = c(0, distance_lon, distance_lon*2)), stringsAsFactors = FALSE, row.names = NULL)
  return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}
#We also need a function to obtain the coordinates of the North arrow:

#
# Result #
#--------#
# Result #
#--------#
# Returns a list containing :
#	- res : coordinates to draw an arrow ;
#	- coordinates of the middle of the arrow (where the "N" will be plotted).
#
# Arguments : #
#-------------#
# scale_bar : result of create_scale_bar() ;
# length : desired length of the arrow ;
# distance : distance between legend rectangles and the bottom of the arrow ;
# dist_units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
create_orientation_arrow <- function(scale_bar, length, distance = 1, dist_units = "km"){
  lon <- scale_bar$rectangle2[1,1]
  lat <- scale_bar$rectangle2[1,2]
  
  # Bottom point of the arrow
  beg_point <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distance, dist.units = dist_units, model = "WGS84")
  lon <- beg_point[1,"long"]
  lat <- beg_point[1,"lat"]
  
  # Let us create the endpoint
  on_top <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = length, dist.units = dist_units, model = "WGS84")
  
  left_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 225, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  right_arrow <- gcDestination(lon = on_top[1,"long"], lat = on_top[1,"lat"], bearing = 135, dist = length/5, dist.units = dist_units, model = "WGS84")
  
  res <- rbind(
    cbind(x = lon, y = lat, xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = left_arrow[1,"long"], y = left_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]),
    cbind(x = right_arrow[1,"long"], y = right_arrow[1,"lat"], xend = on_top[1,"long"], yend = on_top[1,"lat"]))
  
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  
  # Coordinates from which "N" will be plotted
  coords_n <- cbind(x = lon, y = (lat + on_top[1,"lat"])/2)
  
  return(list(res = res, coords_n = coords_n))
}
#The last function enables the user to draw the elements:

#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distance_lon : length of each rectangle ;
# distance_lat : width of each rectangle ;
# distance_legend : distance between rectangles and legend texts ;
# dist_units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec_fill, rec2_fill : filling colour of the rectangles (default to white, and black, resp.);
# rec_colour, rec2_colour : colour of the rectangles (default to black for both);
# legend_colour : legend colour (default to black);
# legend_size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow_length : length of the arrow (default to 500 km) ;
# arrow_distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow_north_size : size of the "N" letter (default to 6).
scale_bar <- function(lon, lat, distance_lon, distance_lat, distance_legend, dist_unit = "km", rec_fill = "white", rec_colour = "black", rec2_fill = "black", rec2_colour = "black", legend_colour = "black", legend_size = 3, orientation = TRUE, arrow_length = 500, arrow_distance = 300, arrow_north_size = 6){
  the_scale_bar <- create_scale_bar(lon = lon, lat = lat, distance_lon = distance_lon, distance_lat = distance_lat, distance_legend = distance_legend, dist_unit = dist_unit)
  # First rectangle
  rectangle1 <- geom_polygon(data = the_scale_bar$rectangle, aes(x = lon, y = lat), fill = rec_fill, colour = rec_colour)
  
  # Second rectangle
  rectangle2 <- geom_polygon(data = the_scale_bar$rectangle2, aes(x = lon, y = lat), fill = rec2_fill, colour = rec2_colour)
  
  # Legend
  scale_bar_legend <- annotate("text", label = paste(the_scale_bar$legend[,"text"], dist_unit, sep=""), x = the_scale_bar$legend[,"long"], y = the_scale_bar$legend[,"lat"], size = legend_size, colour = legend_colour)
  
  res <- list(rectangle1, rectangle2, scale_bar_legend)
  
  if(orientation){# Add an arrow pointing North
    coords_arrow <- create_orientation_arrow(scale_bar = the_scale_bar, length = arrow_length, distance = arrow_distance, dist_unit = dist_unit)
    arrow <- list(geom_segment(data = coords_arrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coords_arrow$coords_n[1,"x"], y = coords_arrow$coords_n[1,"y"], size = arrow_north_size, colour = "black"))
    res <- c(res, arrow)
  }
  return(res)
}


##########################################################################################################

# MAKE SOME PLOTS!!!!

# set the breaks for your color ramp
mybreaks=c(0, 30, 60, 90, 120, 150, 180)


# Make the first map
#tiff(file = "Map_coastline.tiff", width = 36, height = 27, units = "cm", res = 600, compression="lzw") 


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=my_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(Julian day)\n", 
                      breaks = mybreaks) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= my_data, aes(x=longitude, y=latitude, label=code), size=4) +
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 100, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)

#dev.off()


# Make a little map of WA samples

#tiff(file = "Map_WA insert.tiff", width = 36, height = 27, units = "cm", res = 600, compression="lzw") 


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=my_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(days from January 1)\n", 
                      breaks = mybreaks) + 
  labs(element_blank())+
  coord_map(xlim= c(-122, -125),  ylim = c(47,49.1)) +
  geom_text_repel( data= WA_data, aes(x=longitude, y=latitude, label=code), size=4)
  
#dev.off()

# Make a map of BC samples only
#tiff(file = "Map_BC.tiff", width = 36, height = 27, units = "cm", res = 600, compression="lzw") 


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=BC_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(Julian day)\n", 
                      breaks = mybreaks) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= BC_data, aes(x=longitude, y=latitude, label=code), size=4) +
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 100, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)



########################################################
# Make a map with only primary spawners (for powerpoint presentations)

primary_data <- filter(my_data, Spawning != "late spawner", Spawning != "inlet spawner")


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=primary_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(Julian day)\n", 
                      begin = 0, end = 0.60) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= primary_data, aes(x=longitude, y=latitude, label=code), size=4) +
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 100, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)

# Make a map with ONLY late spawners (for powerpoint presentations)

late_data <- filter(my_data, Spawning == "late spawner")


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=late_data, aes(x=longitude, y=latitude, color= days), size = 8, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(Julian day)\n", 
                      begin = 0.7, end = 1) + 
  coord_map(xlim= c(-119, -139),  ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= late_data, aes(x=longitude, y=latitude, label=code), size=4) +
  scale_bar(lon = -136, lat = 47.3, 
            distance_lon = 100, distance_lat = 15, distance_legend = 40, 
            dist_unit = "km", orientation = FALSE)



# Make a map with no late spawners in Salish Sea (for powerpoint presentations)

Salish_primary_data <- filter(my_data, subregion == "Salish Sea",  Spawning != "late spawner" )

ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  geom_point(data=Salish_primary_data, aes(x=longitude, y=latitude, color= days), size = 8, alpha = 0.7) +
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA)) + 
  scale_color_viridis(option="plasma", 
                      name="Date of sampling\n(days from January 1)\n", 
                      begin = 0, end = 0.60) + 
  labs(element_blank())+
  coord_map(xlim= c(-122, -125),  ylim = c(47,49.5)) +
  geom_text_repel( data= Salish_primary_data, aes(x=longitude, y=latitude, label=code), size=4)




# Make the same maps with Different colo schemes!

ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey7", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey27", alpha=0.3)+
  geom_point(data=my_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.9) +
  theme(panel.background = element_rect(fill = "aliceblue")) + 
  scale_color_gradientn(colors = rev(rainbow(4)), name="Date of sampling\n(days from January 1)\n", breaks = mybreaks) + 
  coord_map(xlim= c(-119, -139), ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= my_data, aes(x=longitude, y=latitude, label=code), size=4)


# Different colors


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey7", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey27", alpha=0.3)+
  geom_point(data=my_data, aes(x=longitude, y=latitude, color= days), size = 5, alpha = 0.9) +
  theme(panel.background = element_rect(fill = "aliceblue")) + 
  scale_color_gradientn(colors = (jet(4)), name="Date of sampling\n(days from January 1)\n", breaks = mybreaks) + 
  coord_map(xlim= c(-119, -139), ylim = c(46,60)) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text_repel( data= my_data, aes(x=longitude, y=latitude, label=code), size=3)



#### Random countries


ggplot() +
  geom_polygon(data = USA, aes(x=long, y = lat, group = group), fill=" grey37", alpha=0.3) +
  geom_polygon(data = Canada, aes(x=long, y = lat, group = group), fill=" grey17", alpha=0.3)+
  theme(panel.background = element_rect(fill = "aliceblue"), 
        panel.grid.major = element_line(colour = NA), 
        axis.text=element_text(size=12),
        axis.title =element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) + 
  coord_map(xlim= c(-140, -50),  ylim = c(20,60)) +
  labs(x = "Longitude", y = "Latitude") 


