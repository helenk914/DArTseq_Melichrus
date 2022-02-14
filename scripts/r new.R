#install.packages("rnaturalearthhires", repos = "https://packages.ropensci.org", type = "source")
library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(ggspatial)
library(mapview)
library(rnaturalearth)
library(sp)
library(rnaturalearthhires)
library(rnaturalearth)
library(rnaturalearthdata)


AUS_STATE_shp <- ne_countries(scale = "medium", returnclass = "sf")
class(AUS_STATE_shp)
AUS_STATE_shp <- ne_states(country = 'Australia', returnclass = "sf")
class(AUS_STATE_shp)
ggplot(data = AUS_STATE_shp) +
  geom_sf()



setwd("C:/Users/tshaldoo/OneDrive - University of New England/PhD")


spdat <- read.csv("Baloskion australe.csv")
# remove duplicates
spdat <- unique(spdat)
#make a spatialpoints object
sp <- SpatialPoints(cbind(spdat$longitude, spdat$latitude),
                    proj4string = CRS("+init=epsg:4326"))
# re-add species names
sp$species <- spdat$species
#check the points
mapview(sp) 



#make a new dataset of cities in Australia (google the locations)
AUS_cities <- tribble(
  ~city, ~lat, ~long, 
  "Sydney", -33.8688, 152.2093,
  "QUEENSLAND", -28.5, 145.2093,
  "VICTORIA", -36.8, 145.2093,
  "NEW SOUTH WALES", -32.5000, 145.2631)
#check it looks okay
AUS_cities
#convert those two columns to geometry column with the st_as_sf() function. Google Maps uses the coordinate reference system 4326 (the GPS system).
AUS_cities_geometry <- AUS_cities %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)


AUS_cities1 <- tribble(
  ~city, ~lat, ~long, 
  "Sydney", -33.8688, 152.2093)

AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp, fill="white")+
  coord_sf(xlim = c(140.0, 158.0), ylim = c(-15.5, -38), expand = FALSE)+
  geom_point(data =metamap6, aes(x = longitude, y = latitude), size = 5)+
  geom_point(data =AUS_cities1, aes(x = long, y = lat),
             shape = 1 ,
             size = 5)+
  annotation_scale(location = "br",
                   text_pad = unit(-3, "cm"))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city),
                size = 5, label.size = NA, fill=NA)+
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.8, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle (expression(italic("phebalium glandulosum")))+
  xlab("") +
  ylab("") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 30))

#write.csv(spdat, "phebalium glandulosum.csv")
#unit(-2, "in")

#coord_sf(xlim = c(136.0, 154), ylim = c(-17, -38), expand = FALSE)+

AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp, fill="white")+
  coord_sf(xlim = c(140.0, 158.0), ylim = c(-15.5, -38), expand = FALSE)+
  geom_point(data =metamap6, aes(x = longitude, y = latitude, color = OTU), size = 5)+
  geom_point(data =AUS_cities1, aes(x = long, y = lat),
             shape = 1 ,
             size = 5)+
  annotation_scale(location = "br",
                   text_pad = unit(-3, "cm"))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city),
                size = 5, label.size = NA, fill=NA)+
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.8, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle (expression(italic("Melichrus OTU")))+
  xlab("") +
  ylab("") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 30), panel.background = element_rect(fill = "aliceblue"))

#install.packages("NCmisc")
P <- NCmisc::list.functions.in.file("Short.R")
summary (P)
