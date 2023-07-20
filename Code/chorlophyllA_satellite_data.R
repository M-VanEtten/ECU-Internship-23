library(raster)
library(sp)
library(rerddap)
library(lubridate)
library(readxl)
library(ncdf4)
library(httr)

setwd(dir = "C://Users/fire2/Downloads/ECU Internship '23/")

junk <- 
  GET("https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.htmlTable?chlor_a%5B(2022-12-16T00:00:00Z):1:(2022-12-16T00:00:00Z)%5D%5B(89.97916):1:(-89.97917)%5D%5B(-179.9792):1:(179.9792)%5D"
      , write_disk("cholorA.nc", overwrite = TRUE))


nc = nc_open("cholorA.nc")

