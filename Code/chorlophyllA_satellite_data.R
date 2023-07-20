library(raster)
library(sp)
library(rerddap)
library(lubridate)
library(readxl)
library(ncdf4)
library(httr)
library(dplyr)
library(rerddapXtracto)

setwd(dir = "C://Users/fire2/Downloads/ECU Internship '23/")

# Make sure this URL is the date/time/location range that you want, and that it's the ".nc" version
junk <-
  GET("https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D",
      write_disk("chlorA.nc", overwrite = TRUE))

# Open .nc file
nc = nc_open('chlorA.nc')

# Look at what the variable names are
names(nc$var)

# Extract chlorophyll data
nc_variable1 = nc$var[[1]]
chlor_a = ncvar_get(nc, nc$var[[1]])
dim(chlor_a)

# Extract dates, latitude, and longitude (check these to make sure they look right!)
dates = as.POSIXlt(nc_variable1$dim[[3]]$vals, origin = "1970-01-01", tz="GMT") # I got the time origin by looking at the metadata on ERDDAP ("time_origin")
lon = nc_variable1$dim[[1]]$vals
lat = nc_variable1$dim[[2]]$vals

# ------------------
# Redo

# Import data
prediction_grid <- load(file = "C://KDale/Projects/Phenology/Analysis/PredictionGrids/Anoplopoma fimbria_grid.rdata")
chlor_a = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D"

prediction_grid_subset <- prediction_grid_roms %>% sf::st_drop_geometry()
prediction_grid_subset$dates = paste0(prediction_grid_subset$year, "-", prediction_grid_subset$month, "-", 1, "T00:00:00Z")
prediction_grid_subset$dates = as.POSIXct(prediction_grid_subset$dates)
dates = prediction_grid_subset$dates
lat = prediction_grid_subset$latitude
lon = prediction_grid_subset$longitude

# Redo #2
chlor_a <- griddap(
  datasetx = 'erdMH1chlamday_R2022SQ',
  url = "https://coastwatch.pfeg.noaa.gov/erddap/",
  fields=c('time', 'latitude',  'longitude', 'chlor_a'),
  time = c("2002-07-16T00:00:00Z", "2022-12-16T00:00:00Z"), # This time format needs to be exactly the same as the dataset's
  latitude = c(35,60),
  longitude = c(-134,-117),
  fmt = "nc",
  store = disk()
)

#Create dataframe
chlor.df <-data.frame(longitude=as.numeric(chlor_a$longitude),
                     latitude=as.numeric(chlor_a$latitude),
                     time=strptime(chlor_a$time, "%Y-%m-%dT%H:%M:%S"),
                     chlor_a=as.numeric(chlor_a$chlor_a))

dataInfo <- rerddap::info('erdMH1chlamday_R2022SQ', url = "https://coastwatch.pfeg.noaa.gov/erddap/")
dataInfo

rxtracto(dataInfo, parameter = "chlor_a",
         tcoord = dates,
         xcoord = lon,
         ycoord = lat)
