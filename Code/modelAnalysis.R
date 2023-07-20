# Code written by Martina van Etten

library(readxl)
library(ggplot2)
library(sdmTMB)
library(dplyr)
library(mapping)
library(blockCV)
library(sf)
#head(dataSablefishSDM) #-- to quick view data

#prepare data-------------------------------------------------------------------
setwd(dir = "C://Users/fire2/Downloads/ECU Internship '23/")

dataSablefishSDM <- readxl::read_xlsx("data/sablefish.xlsx")

# create sf object
dataSablefishSF <- st_as_sf(x = dataSablefishSDM, coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>% #geographic reference sys
  st_transform( "EPSG:5070") #set projected coords sys

#download shapefile
NSAmerica <- read_sf("data/North_South_America/North_South_America.shp") %>% st_union() %>%
  st_transform(., crs = "EPSG:5070")



#create mesh -- use make_mesh()-------------------------------------------------
mesh <- make_mesh(dataSablefishSDM, xy_cols = c("X", "Y"), n_knots = 100)

#fit the model -- use sdmTMB()--------------------------------------------------
fitSablefish <- sdmTMB(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled)
                       +as.factor(gearGeneral),
                       data = dataSablefishSDM,
                      # control = sdmTMBcontrol(newton_loops = 1, nlminb_loops = 2),
                       mesh = mesh, family = tweedie(link = "log"),
                       silent = FALSE) #run center of gravity on this model




sanity(fitSablefish)

#create plot--------------------------------------------------------------------

visreg::visreg(fitSablefish, xvar = "sst_scaled", scale ="response")

#plot map of NSAmerica + pos tows
ggplot() +
  geom_sf(data = subset(dataSablefishSF, logN1 > 0), aes(color = logN1)) +
  geom_sf(data = NSAmerica) +
  #super enhances map to view GOA, CAN, NNAMERICA
  xlim(min(dataSablefishSF$X)*1000-1000, max(dataSablefishSF$X)*1000+1000) +
  ylim(min(dataSablefishSF$Y)*1000-1000, max(dataSablefishSF$Y)*1000+1000)

#--(dont worry about this)--model comparison------------------------------------
# cv -- cross validation

# Figure out appropriate block size
sac.ln1 <- cv_spatial_autocor(x = dataSablefishSF, column = "abundance", plot = T)
# plot(sac.ln1$variograms[[1]])
blocksize <- sac.ln1$range

# Several ways of determining folds - spatial, buffer, nearest neighbor
folds.spatial <- cv_spatial(x = dataSablefishSF, size = blocksize,
                            column = "logN1",
                            k = 5,  #number of folds
                            selection = "random",
                            iteration = 50)

dataSablefishSDM$fold = folds.spatial$folds_ids

model1 <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) +s(salinity_scaled) + s(ssh_scaled)
                    + s(distance_from_shore_scaled) +as.factor(gearGeneral),
                    data = dataSablefishSDM,
                    mesh = mesh, family = tweedie(link = "log"),
                    silent = FALSE, fold_ids = "fold" )

# model.SpatioT <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) +s(salinity_scaled) + s(ssh_scaled)
#                     + s(distance_from_shore_scaled) +as.factor(gearGeneral),
#                     data = dataSablefishSDM,
#                     mesh = mesh, family = tweedie(link = "log"),
#                     silent = FALSE, fold_ids = "fold", time = "timeblock",
#                     spatiotemporal = "iid" ) #iid -- identical and independantly destributed (every year separate, not related to one another) )

model.NoSalinity <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled)
                              + s(distance_from_shore_scaled) +as.factor(gearGeneral),
                              data = dataSablefishSDM, spatial = "on",
                              mesh = mesh, family = tweedie(link = "log"),
                              silent = FALSE, fold_ids = "fold" )

#looking for most pos/least neg output
model1$elpd
model.NoSalinity$elpd
model.SpatioT$elpd
#looking for most pos/least neg output
model1$sum_loglik
model.NoSalinity$sum_loglik

#prediction---------------------------------------------------------------------

load(file = "data/Anoplopoma fimbria_grid.rdata")
#try to plot(^^) sf object

pSable <- predict(fitSablefish, newdata = prediction_grid_roms) %>%
  group_by_at(c("region", "year", "latitude", "longitude", "X", "Y"
                ,"timeblock", "sst_roms", "ssh_roms", "month")) %>%
  summarize(est = sum(est), est_rf = mean(est_rf),est_non_rf = mean(est_non_rf)) %>%
  mutate(est_responsescale = fitSablefish$family$linkinv(est))

save(fitSablefish, pSable, file = "Results/prediction_fit_sablefish.rdata")

load(file = "Results/prediction_fit_sablefish.rdata")



#head(pSable)

ggplot(pSable, aes(X, Y, fill = exp(est))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "sqrt")


#sf object based on pridiction data

sablefishPre <- st_as_sf(x = pSable, coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform("EPSG:5070")

ggplot() +
  geom_sf(data = sablefishPre, aes(color = exp(est_responsescale))) +
  geom_sf(data = NSAmerica) +
  facet_wrap(~ timeblock, nrow=1) +
  #super enhances map to view GOA, CAN, NNAMERICA
  xlim(min(dataSablefishSF$X)*1000-1000, max(dataSablefishSF$X)*1000+1000) +
  ylim(min(dataSablefishSF$Y)*1000-1000, max(dataSablefishSF$Y)*1000+1000)










