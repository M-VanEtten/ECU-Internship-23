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

# Import prediction grid
load(file = "C://KDale/Projects/Phenology/Analysis/PredictionGrids/Anoplopoma fimbria_grid.rdata")

# Set URL for the satellite data file we're interested in
chlor_a = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D"

# Subset to date range covered in satellite data
prediction_grid = sf::st_drop_geometry(prediction_grid_roms) %>%
  mutate(dates = paste0(year, "-", month, "-", 1, "T00:00:00Z")) %>%
  mutate(dates = as.POSIXct(dates)) %>%
  subset(dates > "2002-07-16T00:00:00Z" & dates < "2022-12-16T00:00:00Z")

# Access satellite data
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

dataInfo <- rerddap::info('erdMH1chlamday_R2022SQ', url = "https://coastwatch.pfeg.noaa.gov/erddap/")
dataInfo

extract <- rxtracto(dataInfo, parameter = "chlor_a",
         tcoord = prediction_grid$dates,
         xcoord = prediction_grid$longitude,
         ycoord = prediction_grid$latitude)

save(extract, file = "chlor_a.rdata")
load("chlor_a.rdata")

prediction_grid$chlor_a = extract$`mean chlor_a`
save(prediction_grid, file = "Anoplopoma fimbria_grid.rdata")

