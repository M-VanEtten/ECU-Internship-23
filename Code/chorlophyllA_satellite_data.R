library(raster)
library(sp)
library(rerddap)
library(lubridate)
library(readxl)
library(ncdf4)
library(httr)
library(dplyr)

setwd(dir = "C://Users/fire2/Downloads/ECU Internship '23/")

# Make sure this URL is the date/time/location range that you want, and that it's the ".nc" version
junk <-
  GET("https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D",
      write_disk("chlorA.nc", overwrite = TRUE))

# Open .nc file
nc = nc_open("chlorA.nc")

# Look at what the variable names are
names(nc$var)

# Extract chlorophyll data
chlor_a = ncvar_get(nc, nc$var[[1]])



# ------------------
# Redo

# Import data
prediction_grid <- load(file = "C://KDale/Projects/Phenology/Analysis/PredictionGrids/Anoplopoma fimbria_grid.rdata")
chlor_a = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D"

dates <- prediction_grid_roms %>% sf::st_drop_geometry() %>% subset(year >= 2002 & month >= 8)
dates <- paste(dates$year, dates$month, 1, sep = "-")
dates <- format(as.Date(dates), "%Y-%m-%d")
lat = prediction_grid_roms$latitude
lon = prediction_grid_roms$longitude

tot = rep(NA, 4)
pb <- txtProgressBar(min = 0, max = nrow(sablefish), char = "=", style = 3)

for (i in 1:nrow(sablefish)) {

  #this is where the URL is built:
  url = paste(chlor_a, "[(", dates[i], "):1:(", dates[i], ")][(", lat[i], "):1:(", lat[i], ")][(", lon[i], "):1:(", lon[i], ")]", sep = "")

  new = read.csv(url, skip = 2, header = FALSE)
  tot = bind_rows(tot, new)

  setTxtProgressBar(pb, i)

}
close(pb) # Close progress bar
