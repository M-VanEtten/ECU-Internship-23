library(raster)
library(sp)
library(rerddap)
library(lubridate)
library(readxl)
library(ncdf4)
library(httr)
library(dplyr)
library(rerddapXtracto)

linkChlorophyll <- function() {
# Import prediction grid
load(file = "data/Anoplopoma fimbria_grid.rdata")

# Set URL for the satellite data file we're interested in
chlor_a = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D"

# Subset to date range covered in satellite data
prediction_grid_chlorA = sf::st_drop_geometry(prediction_grid_roms) %>%
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
  store = disk())

# Get info about the satellite dataset
dataInfo <- rerddap::info('erdMH1chlamday_R2022SQ', url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# Extract data at our locations (this takes a few minutes -- can reload using lines below)
extract <- rxtracto(dataInfo, parameter = "chlor_a",
         tcoord = prediction_grid_chlorA$dates,
         xcoord = prediction_grid_chlorA$longitude,
         ycoord = prediction_grid_chlorA$latitude)

# Save/load extraction file
save(extract, file = "chlor_a.rdata")
load("chlor_a.rdata")

# Add chlorophyll data to our prediction grid
prediction_grid_chlorA$chlor_a = extract$`mean chlor_a`
prediction_grid_chlorA = subset(prediction_grid_chlorA, chlor_a < 10) # remove outliers

# Save grid
save(prediction_grid_chlorA, file = "data/Anoplopoma fimbria_grid_chlorA.rdata")
}
