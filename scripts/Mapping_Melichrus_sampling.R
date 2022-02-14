################################################################################
## Mapping Melichrus sampling regime
################################################################################

## Install packages ##
# install.packages("rnaturalearthhires", repos = "https://packages.ropensci.org", type = "source") 
# install.packages("rnaturalearthdata")

## Load packages ##

library(dplyr)
library(ggplot2)
library(sf)
library(raster)
library(ggspatial)
library(mapview)
library(sp)
library(rnaturalearthhires)
library(rnaturalearth)
library(rnaturalearthdata)

## Download base map of Australia ##
AUS_STATE_shp <- ne_states(country = 'Australia', returnclass = "sf")
class(AUS_STATE_shp)
ggplot(data = AUS_STATE_shp) + geom_sf() 

## Read in metadata ##

metamap0 <- read.csv("data/metadata.csv", stringsAsFactors = F)

## Remove ARC870 ##

metamap <- metamap0[-c(551),]

## Make new dataframe with only coordinates, OTU, collection number, and collector ID ##

metamap1 <- dplyr::select(metamap, decimal_longitude, decimal_latitude, scientific_name_OTU, voucher_herbarium_collector_id, voucher_collection_number)
colnames(metamap1) <- c("longitude", "latitude", "OTU","collectorID", "collectionNumber")


## Remove population duplicates ##
metamap2 <- unique(metamap1)

## Create a spatial points object ##
MELsp <- SpatialPoints(cbind(metamap2$longitude, metamap2$latitude),
                       proj4string = CRS("+init=epsg:4326"))

## Add species names (OTU) to spatial points objects

MELsp$OTU <- metamap2$OTU

## check the points ##

mapview(MELsp) 

## Make a dataframe with the names of States and landmark cities and location for label ##

AUS_cities <- tribble(
  ~city, ~lat, ~long, 
  "Sydney", -33.8688, 152.2093,
  "QUEENSLAND", -22.5, 143.5,
  "VICTORIA", -36.8, 142.6,
  "NEW SOUTH WALES", -30, 144.5)

# Convert the lat and long columns to geometry column with the st_as_sf() function.
# Google Maps uses the coordinate reference system 4326 (the GPS system).
AUS_cities_geometry <- AUS_cities %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

## Create a separate file to produce circle at location of Sydney ##
AUS_cities1 <- tribble(
  ~city, ~lat, ~long, 
  "Sydney", -33.8688, 151.1)

## Format map ## 
AUS_cities_geometry
ggplot() +
  geom_sf(data=AUS_STATE_shp, fill="white")+
  coord_sf(xlim = c(141.0, 158.2), ylim = c(-15.5, -38), expand = FALSE)+
  geom_point(data =metamap2, aes(x = longitude, y = latitude, color = OTU), size = 5)+
  geom_point(data =AUS_cities1, aes(x = long, y = lat),
             shape = 1 ,
             size = 2)+
  annotation_scale(location = "br",
                   text_pad = unit(0.1, "cm"))+
  geom_sf_label(data =AUS_cities_geometry, aes(label = city),
                size = 3, label.size = NA, fill=NA)+
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.2, "in"), 
                         pad_y = unit(0.3, "in"), 
                         style = north_arrow_fancy_orienteering) +
  ggtitle (expression(italic("Melichrus")))+xlab("") +
  ylab("") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 30), panel.background = element_rect(fill = "aliceblue"))




## To do ##
# Create palette
# Apply palette
# Order specimens so the palette is applied correctly and they appear in order that you want in legend
# Fade State Labels













