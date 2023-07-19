#install packages---------------------------------------------------------------
install.packages("ncdf4")
install.packages("httr")

library(ncdf4)
library(httr)

#run the URL--------------------------------------------------------------------

junk <- 
  GET("https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v1_0_monthly.nc?analysed_sst[(2018-01-01T12:00:00Z):1:(2018-12-01T12:00:00Z)][(17):1:(30)][(195):1:(210)]"
      , write_disk("sst.nc", overwrite = TRUE))

# importing the downloaded data-------------------------------------------------

  #-- opening file
nc=nc_open("sst.nc")

  #--examine variables included within the dataset
names(nc$var)

 #--extract analysed_sst
v1=nc$var[[1]]
sst=ncvar_get(nc,v1)

 #--examine the structure of sst
dim(sst) #dataset is a 3-D array w/ 301 rows(longitude), 261 columns(latitude), 12(time steps)

 #--get the dates for each time step
dates=as.POSIXlt(v1$dim[[3]]$vals, origin = "1970-01-01", tz = "GMT")

 #--get longitude and latitude values
lon = v1$dim[[1]]$vals
lat = v1$dim[[2]]$vals

 #--close netcdf file adn remove unneeded data and files
nc_close(nc)
rm(junk,v1)
file.remove("stt.nc")

#working with the extracted data------------------------------------------------

#--creating a map for one time step
 #--set some color breaks
h = hist(sst[,,1], 100, plot = FALSE)
breaks = h$breaks
n = length(breaks)-1

 #--define a color palette
jet.colors <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F",
                                 "yellow", "#FF7F00", "red", "#7F0000"))

 #set color scale using jet.color palette
c = jet.colors(n)

 #prepare graphic window: left side for map, right side for color scale
layout(matrix(c(1,2,3,0,4,0), nrow = 1, ncol = 2), widths = c(5,1), 
       heights = 4)
layout.show(2)
par(mar = c(3,3,3,1))

 #--plot the SST map
image(lon, lat, sst[,,1], col = c, breaks = breaks, xlabs = "X", ylabs = "Y", 
      axes = TRUE, xaxs = "i", yaxs = "i", asp = 1, main = paste("Monthly SST",
                                                                 dates[1]))
 #--example on how to add points
points(202:205, rep(26,4), pch = 20, cex = 2)

 #example on how to add contour (condidered a new plot, not a feature, need to 
 #                                 overlay on top of SST map)
par(new = TRUE)
contour(lon, lat, sst[,,1], levels = 20, xaxs = "i", yaxs = "i", labcex = 0.8, 
        vfont = c("sans serif", "bold", axes = FALSE, asp = 1))
 
 #--plot color scale --  
 #-- paste this into the working directory(vv)
# image.scale <- function(z, zlim, col = heat.colors(12),
#                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
#   if(!missing(breaks)){
#     if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
#   }
#   if(missing(breaks) & !missing(zlim)){
#     breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
#   }
#   if(missing(breaks) & missing(zlim)){
#     zlim <- range(z, na.rm=TRUE)
#     zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
#     zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
#     breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
#   }
#   poly <- vector(mode="list", length(col))
#   for(i in seq(poly)){
#     poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
#   }
#   xaxt <- ifelse(horiz, "s", "n")
#   yaxt <- ifelse(horiz, "n", "s")
#   if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
#   if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
#   if(missing(xlim)) xlim=XLIM
#   if(missing(ylim)) ylim=YLIM
#   plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
#   for(i in seq(poly)){
#     if(horiz){
#       polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
#     }
#     if(!horiz){
#       polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
#     }
#   }

 #to make the acutal color scale
par(mar = c(3,1,3,3))
source("scale.R")
image.scale(sst[..1], col = c, breaks = breaks, horiz = FALSE, yaxt = "n", 
            xlabs = "X", , ylabs = "Y", main = "SST")
axis(4, las = 1) 
box()

#plotting a time series---------------------------------------------------------

I= which(lon>=200 & lon<=206)
J = which(lat>= 24 & lat <= 26) 
sst2 = sst[I,J,]

n = dim(sst2) [3]

res = rep(NA, n) 
for(i in 1:n)
  res[i] = mean(sst2[,,i], na.rm = TRUE)
plot(1:n, res, axes = FALSE, type = "o", pch = 20, xlab = "Month", 
     ylab = "SST (C)")
axis(2)
axis(1,1:n, format(dates, "%m"))
box()


#creating a map of average SST over a year--------------------------------------

sst.yr = apply(sst[,,1:12], c(1,2), mean, na.rm = TRUE)

h = hist(sst.yr, 100, plot = FALSE)
breaks = h$breaks
n = length(breaks) -1
c = jet.colors(n)


layout(matrix(c(1,2,3,0,4,0), nrow = 1, ncol = 2), widths = c(5,1), 
              heights = 4)
layout.show(2)

par(mar = c(3,3,3,1))

image(lon, lat, sst.yr, col = c, breaks = breaks, xlab = "X", ylabs = "Y", 
      axes = TRUE, xaxs = "i", yaxs = "i", asp = 1, main = paste("Mean SST",
                                                                 format(dates[1], 
                                                                        "%Y/%m/%d")))

par(mar = c(3,1,3,3))
image.scale(sst.yr, col = c, breaks = breaks, horiz = FALSE, yaxt = "n", 
            xlabs = "x", ylabs = "y", main = "SST")
axis(4)
box()










