# Code written by Martina van Etten for GSK internship

library(readxl)
library(openxlsx)
library(ggplot2)
library(sdmTMB)
library(dplyr)
library(mapping)
library(blockCV)
library(sf)
library(tidyr)
#head(dataSablefishSDM) #-- to quick view data

#prepare data (function to combine Bering Sea and CCE -------------------------
createSpeciesDataset <- function() {

  species = "Anoplopoma fimbria"

  all_tows_roms <- read.csv("C://KDale/Projects/Phenology/Data/AllTows_200nm.csv")
  beringSeaCTD <- read.csv("data/BASIS_CTD_48e9_1d95_4312.csv") %>%
    group_by_at(c("station_id", "latitude", "longitude", "time")) %>%
    summarize(surface_temp_oC = mean(primary_temperature), surface_sal_psu = mean(primary_salinity), chlorophyll_a = mean(chlorophyll_a))

  beringSea <- read.csv("data/BASIS_FishCatch_a759_d3cb_ab11.csv") %>%
    subset(., gearGeneral == "CobbMWT") %>%
    select(-c(latitude, longitude, time)) %>%
    merge(., beringSeaCTD, all.y = T, by = c("station_id")) %>%
    rename(larvae_count = total_catch_num) %>%
    tidyr::separate(., col = time, into = c("date", "time"), sep = c("T")) %>%
    mutate(latitude = coalesce(eq_latitude, latitude), longitude = coalesce(eq_longitude, longitude)) %>%
    mutate(gearGeneral = "CobbMWT", scientific_name = "Anoplopoma fimbria", month = lubridate::month(date), year = lubridate::year(date),  day = lubridate::day(date), program = "BASIS") %>%
    select(., c("program", "latitude", "longitude", "year", "month", "day", "time", "gearGeneral", "larvae_count", "surface_temp_oC", "surface_sal_psu"))

  speciesData <- read.csv(file = "c://KDale/Projects/Phenology/Data/AllCruises_Combined_200nm.csv") %>%
    subset(scientific_name == "Anoplopoma fimbria" & year >= 1995 & year <= 2019) %>%
    merge(., all_tows_roms, all.y = TRUE) %>%
    bind_rows(., beringSea) %>%
    subset(gearGeneral != "MOCNESS") %>%
    mutate(., presence = 1) %>%
    mutate(., abundance_scaled = coalesce(
      scale(larvae_10m2, center = F)[,1],
      scale(larvae_m3, center = F)[,1],
      scale(larvae_1000m3, center = F)[,1],
      scale(larvae_count, center = F)[,1]
    )) %>%
    mutate(., abundance = coalesce(
      larvae_10m2,
      larvae_m3,
      larvae_1000m3,
      larvae_count
    ))

  # Liknk ROMS data.
  source("Code/linkRoms.R")
  speciesData <- linkroms(speciesData)

  # Get maximum extents of positive tows
  maxLat = max(speciesData$latitude)
  minLat = min(speciesData$latitude)
  maxLon = max(speciesData$longitude)
  minLon = min(speciesData$longitude)

  # Replace all NA abundance values with zero (these are true zeroes)
  speciesData <-  mutate(speciesData, presence = replace_na(presence, 0),
                         abundance = replace_na(abundance, 0),
                         abundance_scaled = replace_na(abundance_scaled, 0),
                         scientific_name = species) %>%
    mutate(sst = coalesce(surface_temp_oC, sst_roms), salinity = coalesce(surface_sal_psu, salinity_roms)) %>%
    subset(sst != 0 & salinity != 0 & sst <= 40 & salinity <= 40) %>% # remove lower and upper outliers
    subset(!is.na(sst) & !is.na(ssh_roms) & !is.na(salinity)) %>% # remove any tows with missing info
    mutate(., sst_scaled = scale(sst)[,1], ssh_scaled = scale(ssh_roms)[,1], salinity_scaled = scale(salinity)[,1]) %>%  # Center and scale enviro data
    mutate(., logN1 = log(abundance_scaled+1))

  # Add Conus Albers coordinates
  coords= cbind(speciesData$longitude, speciesData$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  speciesData$X = albert_equal_area[,1]
  speciesData$Y = albert_equal_area[,2]

  # Add time block
  timeblocks <- read_xlsx("C://KDale/Projects/Phenology/Data/timeblocks.xlsx", sheet =  1)
  speciesData <- merge(speciesData, timeblocks, by = "year")
  speciesData$timeblock = factor(speciesData$timeblock, levels = c("1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019"))

  openxlsx::write.xlsx(x = speciesData, file = "data/sablefish.xlsx")

  return(speciesData)
}
#------------------------------------------------------------------------------
# Run function above (only necessary to do once unless you have new data)
#dataSablefishSDM <- createSpeciesDataset()

# Load data
dataSablefishSDM <- read_xlsx("data/sablefish.xlsx", sheet = 1) %>%
  dplyr::select(!timeblock)
timeblocks <- read_xlsx("C://KDale/Projects/Phenology/Data/timeblocks.xlsx", sheet =  1)
dataSablefishSDM <-subset(dataSablefishSDM, latitude > 35 & year <= 2018 & year >= 1995) %>%
  merge(., timeblocks, by = "year") %>%
  mutate(timeblock_factor = as.numeric(as.factor(timeblock)), timeblock = as.factor(timeblock))

openxlsx::write.xlsx(dataSablefishSDM, "data/sablefish.xlsx")

# Subset data to higher latitudes

#create mesh -- use make_mesh()-------------------------------------------------
mesh <- make_mesh(dataSablefishSDM, xy_cols = c("X", "Y"), n_knots = 200)

#fit the model -- use sdmTMB()--------------------------------------------------
fitSablefish <- sdmTMB(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled) + s(salinity_scaled) + as.factor(gearGeneral) ,
                       spatial = "on",
                       data = dataSablefishSDM,
                       # control = sdmTMBcontrol(newton_loops = 1, nlminb_loops = 2),
                       mesh = mesh, family = tweedie(link = "log"), spatiotemporal = "iid", time = "timeblock",
                       silent = FALSE) #run center of gravity on this model

sanity(fitSablefish)

visreg::visreg(fitSablefish, xvar = "sst_scaled")
visreg::visreg(fitSablefish, xvar = "ssh_scaled")
visreg::visreg(fitSablefish, xvar = "salinity_scaled")

#create plot--------------------------------------------------------------------

#create sf object
dataSablefishSF <- st_as_sf(x = dataSablefishSDM, coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>% #geographic reference sys
  st_transform( "EPSG:5070") #set projected coords sys

#download shapefile
NSAmerica <- read_sf("data/North_South_America/North_South_America.shp") %>% st_union() %>%
  st_transform(., crs = "EPSG:5070")

#plot map of NSAmerica + pos tows
jpeg(filename = "Results/sablefish_positive_tows.jpg", width = 4, height = 6, units = "in", res = 400)
ggplot() +
  geom_sf(data = subset(dataSablefishSF, logN1 > 0), aes(color = logN1)) +
  geom_sf(data = NSAmerica, fill = "gray70") +
  theme_classic() +
  #crops map to view GOA, CAN, NSAMERICA
  xlim(min(dataSablefishSF$X)*1000-1000, max(dataSablefishSF$X)*1000+1000) +
  ylim(min(dataSablefishSF$Y)*1000-1000, max(dataSablefishSF$Y)*1000+1000)
dev.off()

# #--(dont worry about this)--model comparison------------------------------------
# See the model comparison code!
# # cv -- cross validation
#
# # Figure out appropriate block size
# sac.ln1 <- cv_spatial_autocor(x = dataSablefishSF, column = "abundance", plot = T)
# # plot(sac.ln1$variograms[[1]])
# blocksize <- sac.ln1$range
#
# # Several ways of determining folds - spatial, buffer, nearest neighbor
# folds.spatial <- cv_spatial(x = dataSablefishSF, size = blocksize,
#                             column = "logN1",
#                             k = 5,  #number of folds
#                             selection = "random",
#                             iteration = 50)
#
# dataSablefishSDM$fold = folds.spatial$folds_ids
#
# model1 <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) +s(salinity_scaled) + s(ssh_scaled)
#                     + s(distance_from_shore_scaled) +as.factor(gearGeneral),
#                     data = dataSablefishSDM,
#                     mesh = mesh, family = tweedie(link = "log"),
#                     silent = FALSE, fold_ids = "fold" )
#
# # model.SpatioT <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) +s(salinity_scaled) + s(ssh_scaled)
# #                     + s(distance_from_shore_scaled) +as.factor(gearGeneral),
# #                     data = dataSablefishSDM,
# #                     mesh = mesh, family = tweedie(link = "log"),
# #                     silent = FALSE, fold_ids = "fold", time = "timeblock",
# #                     spatiotemporal = "iid" ) #iid -- identical and independantly destributed (every year separate, not related to one another) )
#
# model.NoSalinity <- sdmTMB_cv(formula = logN1 ~ s(sst_scaled) + s(ssh_scaled)
#                               + s(distance_from_shore_scaled) +as.factor(gearGeneral),
#                               data = dataSablefishSDM, spatial = "on",
#                               mesh = mesh, family = tweedie(link = "log"),
#                               silent = FALSE, fold_ids = "fold" )
#
# #looking for most pos/least neg output
# model1$elpd
# model.NoSalinity$elpd
# model.SpatioT$elpd
# #looking for most pos/least neg output
# model1$sum_loglik
# model.NoSalinity$sum_loglik

#prediction---------------------------------------------------------------------
source("Code/createPredictionGrid.R")
# Create prediction grid -- this takes a long time! It will save the grid as part of the function
#createPredictionGrid(data = dataSablefishSDM, species = "Anoplopoma fimbria", path = "data/Anoplopoma fimbria_grid.rdata")
load(file = "data/Anoplopoma fimbria_grid.rdata")

# Subset prediction grid to latitudes >35
prediction_grid_roms = subset(prediction_grid_roms, latitude > 35)

prediction_grid_roms <- prediction_grid_roms %>% group_by_at(c("year", "X", "Y", "region", "latitude", "longitude", "timeblock", "timeblock_factor", "gearGeneral")) %>%
  summarize(sst_scaled = mean(sst_scaled), sst_roms = mean(sst_roms), ssh_roms = mean(ssh_roms), salinity_roms = mean(salinity_roms), ssh_scaled = mean(ssh_scaled), salinity_scaled = mean(salinity_scaled)) %>%
  mutate(timeblock_factor = as.numeric(as.factor(timeblock)))

prediction_grid.sf = st_as_sf(prediction_grid_roms, coords = c("longitude", "latitude"))

# Predict for the current time period on a full grid, then summarize
pSable <- predict(fitSablefish, newdata = prediction_grid_roms) %>%
  group_by_at(c("region", "year", "latitude", "longitude", "X", "Y", "timeblock")) %>%
  summarize(est = mean(est), est_rf = mean(est_rf),est_non_rf = mean(est_non_rf), mean_sst = mean(sst_roms), mean_ssh = mean(ssh_roms), mean_salinity = mean(salinity_roms)) %>%
  ungroup()

# Chlorophyll -- need to revise these lines/get satellite data
#mutate(chlor_a_scaled = scale(mean_chlor_a)[,1], est_scaled = scale(exp(fitSablefish$family$linkinv(est)-1))[,1]) %>% # center and scale chlorophyll
# mutate(est_chlor_diff = est_scaled - chlor_a_scaled) # calculate difference between predicted abundance and chlor_a

#sf object based on prediction data
pSable.sf <- st_as_sf(x = pSable, coords = c("longitude", "latitude")) %>%
  st_set_crs(4326) %>%
  st_transform("EPSG:5070")

# Create prediction object (note the "return_tmb_object = T")
pSable.obj <- predict(fitSablefish, newdata = prediction_grid_roms, return_tmb_object = T)

# Write/read sablefish prediction results -------------------------------------
save(fitSablefish, pSable, pSable.obj, pSable.sf, file = "Results/prediction_objects_sablefish.rdata")
load(file = "Results/prediction_objects_sablefish.rdata")

# FIGURES ----------------------------------------------------------------------
# Plot biomass per timeblock on the original abundance scale
ggplot() +
  geom_sf(data = pSable.sf, aes(color = exp(fitSablefish$family$linkinv(est)-1))) +
  geom_sf(data = NSAmerica, fill = "gray70", color = "black") +
  facet_wrap(~ timeblock, nrow=1) +
  scale_color_gradient("Estimate", low = "khaki1", high = "orchid4") +
  theme_classic() +
  #super enhances map to view GOA, CAN, NNAMERICA
  xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
  ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000)

# Plot match-mismatch difference (hint: use code above and change "color = ...")
# ggplot() +
#   geom_sf(data = pSable.sf, aes(color = est_chlor_diff)) +
#   geom_sf(data = NSAmerica) +
#   facet_wrap(~ timeblock, nrow=1) +
#   theme_classic() +
#   scale_color_gradient2("Mismatch", low = "darkred", mid = "green3", high = "darkred")+
#   #super enhances map to view GOA, CAN, NNAMERICA
#   xlim(min(prediction_grid_roms$X)*1000-1000, max(prediction_grid_roms$X)*1000+1000) +
#   ylim(min(prediction_grid_roms$Y)*1000-1000, max(prediction_grid_roms$Y)*1000+1000)

# CENTER OF GRAVITY ------------------------------------------------------------

# Run center of gravity function provided by sdmTMB
cog <- get_cog(pSable.obj, format = "long")

# Get this into a format we can do things with
cog.wide <-  pivot_wider(cog,  names_from = "coord", values_from = c("est", "lwr", "upr", "se"))

# Plot
ggplot(cog.wide, aes(est_X, est_Y, colour = timeblock)) +
  geom_pointrange(aes(xmin = lwr_X, xmax = upr_X)) +
  geom_pointrange(aes(ymin = lwr_Y, ymax = upr_Y)) +
  theme_classic() +
  scale_color_brewer(palette = "RdPu")

# Calculate center of gravity manually ----
calculateCenterOfGravity <- function(p, species) {

  # Center of gravity per year in relation to environmental variables
  years = unique(p$year)
  cog.year = data.frame(matrix(ncol = 0, nrow = length(years)), year = NA, timeblock = NA, cog.lat = NA, cog.lon = NA, upr_y = NA, lwr_y = NA, upr_x = NA, lwr_x = NA,  avgsst = NA, avgssh = NA, avgsal= NA)

  for (i in 1:length(years)) {
    yearSubset = subset(p, year == years[i])
    cog.year$cog.lat[i] = sum(yearSubset$est*yearSubset$latitude)/sum(yearSubset$est)
    cog.year$cog.lon[i] = sum(yearSubset$est*yearSubset$longitude)/sum(yearSubset$est)
    cog.year$upr_y[i] = max(yearSubset$latitude)
    cog.year$lwr_y[i] = min(yearSubset$latitude)
    cog.year$upr_x[i] = max(yearSubset$longitude)
    cog.year$lwr_x[i] = min(yearSubset$longitude)
    cog.year$avgsst[i] = mean(yearSubset$mean_sst)
    cog.year$avgssh[i] = mean(yearSubset$mean_ssh)
    cog.year$avgsal[i] = mean(yearSubset$mean_salinity)
    cog.year$year[i] = years[i]
    cog.year$species[i] = species
    cog.year$timeblock[i] = as.character(yearSubset$timeblock[1])
  }

  return(cog.year)

}

# Calculate center of gravity using our custom function above
cog <- calculateCenterOfGravity(pSable, species = "Anoplopoma fimbria") %>%
  mutate(cog = as.character(timeblock))

# Plot each year, color-coded by timeblock
ggplot(cog, aes(cog.lon, y = cog.lat, color = timeblock)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "RdPu") +
  theme_classic(base_size = 14) +
  labs(x = "Longitude", y = "Latitude")

