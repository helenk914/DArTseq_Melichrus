##############################################################
################# Species Map ################################
##############################################################

#install spatial packages
install.packages("plyr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("rgdal")
install.packages("tmap")
install.packages("ggmap")
install.packages("sf")
install.packages("broom")
install.packages("tidyverse")
install.packages("readxl")
install.packages("raustats")
install.packages("purrr")
install.packages("Census2016")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("ggspatial")
install.packages("rJava")
install.packages("dismo")
install.packages("rjasonlite")
install.packages("mapdata")
install.packages("rgbif")
install.packages("mapview")
install.packages("scrubr")
install.packages("sp")



#load spatial packages
library(plyr)
library(dplyr)
library(ggplot2)
library(rgdal)
library(tmap)
library(ggmap)
library(sf)
library(ggspatial)
library(rlang)
library(broom)
library(tidyverse)
library(readxl)
library(raustats)
library(purrr)
library("Census2016")
library(rnaturalearth)
library(rnaturalearthdata)
library(dismo)
library(jsonlite)
library(mapdata)
library(raster)
library(mapview)
library(rgbif)
library(maptools)
library(rgbif)
library(scrubr)
library(sp)





#set the working directory 
getwd()



#################### Download maps in ESRI Shapefile format from here ##########################
## https://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/1270.0.55.001July%202016?OpenDocument ##



#####################GET the data of your spices if you do not have Data########################
# Pull records on a family or ... from GBIF
key <- name_suggest(q='Phebalium glandulosum', rank='species')$data$key[1]
# Show metadata of records for phebalium in the database
occ_search(taxonKey=key, limit=0)$meta$count
# Pulls your data from GBIF, limit to 200 records as an example dataset
spdat <- occ_search(taxonKey = key, return = "data", limit = 963)
#pull out the data file
spdat <- spdat$data
#view the data that was returned
View(spdat)

#write.csv(spdat, "Phebalium glandulosum.csv")

###############################If you have the Data#################################
list.files()
dat <- read.csv("Phebalium glandulosum.csv")
# keep only coordinates and species name
#spdat <- select(dat, coords.x1, coords.x2, species)
#colnames(spdat) <- c("longitude", "latitude", "species")

spdat <- dplyr::select(dat, decimalLongitude, decimalLatitude, scientificName, catalogNumber)
colnames(spdat) <- c("longitude", "latitude", "species","catalogNumber")

# remove NA values
spdat <- na.omit(spdat)

###clean the coordinates##
#dat<-dframe(spdat) %>% coord_impossible()
#dat<-dframe(spdat) %>% coord_unlikely()


# remove duplicates
spdat <- unique(spdat)
#make a spatialpoints object
sp <- SpatialPoints(cbind(spdat$longitude, spdat$latitude),
                    proj4string = CRS("+init=epsg:4326"))
# re-add species names
sp$species <- spdat$species
#check the points
mapview(sp) 


############################# first method ################################





#dat.all <- read.csv("glandulosum.csv") #copy/paste the data in the function

#spdat <- select(dat.all, decimalLongitude, decimalLatitude, scientificName, catalogNumber)
#colnames(spdat) <- c("longitude", "latitude", "species","catalogNumber")

#spdat <- na.omit(spdat)

# clean the coordinates
#dat<-dframe(spdat) %>% coord_impossible()
#dat<-dframe(spdat) %>% coord_unlikely()

#make a spatialpoints object
#sp <- SpatialPoints(cbind(spdat$longitude, spdat$latitude),
                    proj4string = CRS("+init=epsg:4326"))

#sp$species <- spdat$species

# quickly map and check our points to make sure there is 
# nothing out of the ordinary (points in the ocean/arctic etc.)
#mapview(sp) 



######################### Second method #######################

#import a shapefile of state boundaries
AUS_STATE_shp <- read_sf("C:/Users/tshaldoo/Documents/R example/ASGS","STE_2016_AUST")

#list.files() # find the file we just made called "dataframe_file.csv"   
#dat <- read.csv("Phebalium glandulosum subsp. angustifolium Paul G.Wilson.csv") #copy/paste the data in the function
#View(dat) #

# keep only coordinates and species name
#spdat <- select(dat, coords.x1, coords.x2, species)
#colnames(spdat) <- c("long", "lat", "species")
# remove NA values
#spdat <- na.omit(spdat)

#make a new dataset of cities in Australia (google the locations)
AUS_cities <- tribble(
  ~city, ~lat, ~long, 
  "Brisbane",-27.4698, 153.0251,
  "Sydney", -33.8688, 151.2093,
  "Melbourne", -37.8136, 144.9631)
#check it looks okay
AUS_cities
#convert those two columns to geometry column with the st_as_sf() function. Google Maps uses the coordinate reference system 4326 (the GPS system).
AUS_cities_geometry <- AUS_cities %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

#check it looks right
AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp)+
  geom_point(data =spdat, aes(x = longitude, y = latitude))+
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle("Phebalium glandulosum subsp. angustifolium Paul G.Wilson") +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#check it looks right with cities
AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp)+
  geom_point(data =spdat, aes(x = longitude, y = latitude))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city))+
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle("Phebalium glandulosum subsp. angustifolium Paul G.Wilson") +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()


#################NSW################
AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp)+
  coord_sf(xlim = c(140, 155), ylim = c(-25, -40), expand = FALSE)+
  geom_point(data =spdat, aes(x = longitude, y = latitude, col='red'))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city))+
  annotation_north_arrow(location = "tl", 
                         which_north = "true", 
                         pad_x = unit(0.1, "in"), 
                         pad_y = unit(0.1, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle("Phebalium glandulosum subsp. angustifolium") +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()



####################more than sp################   
AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp)+
  coord_sf(xlim = c(140, 155), ylim = c(-25, -40), expand = FALSE)+
  geom_point(data =spdat, aes(x = longitude, y = latitude, col='red'))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city))+
  annotation_north_arrow(location = "tl", 
                         which_north = "true", 
                         pad_x = unit(0.1, "in"), 
                         pad_y = unit(0.1, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle("Phebalium glandulosum subsp. angustifolium") +
  xlab("Longitude") +
  ylab("Latitude") + 
  scale_color_discrete(name = "sp", labels = c("sp1", 
                                                 "sp2",
                                                     "sp3"))+
theme_bw()










####https://data.gov.au/dataset/ds-dga-a1b278b1-59ef-4dea-8468-50eb09967f18/details



#import a shapefile of state boundaries 
#AUS_STATE_shp <- read_sf("C:/Users/tshaldoo/Documents/R example/ASGS","MB_2016_NSW")
#AUS_STATE_shp <- read_sf("C:/Users/tshaldoo/Documents/R example/ASGS","NSW_STATE_POLYGON_shp_GDA2020")

#AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp)+
  geom_point(data =spdat, aes(x = longitude, y = latitude, col = "green"))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city))+
  annotation_north_arrow(location = "tr", 
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle("Phebalium glandulosum subsp. angustifolium Paul G.Wilson") +
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#plot(AUS_STATE_shp)
