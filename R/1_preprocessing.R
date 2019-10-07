rm(list = ls())
#required packages
library("data.table")
library('rgdal')
library("rgeos")
library("matrixStats")
require("rgdal")
require("sp")
library('raster')
require(WorldClimTiles)

source(file.path(".", "R", "1_functions.R"))
country_abbr <- "aus" #this script needs to be run for aus, til and vnm
layer_path <- file.path(".", "RData")

#To be done once initially
mask <- getData("GADM", country = country_abbr, level = 0, path = layer_path)
mask <- gSimplify(mask, tol = 0.00833)

if(country_abbr == "vnm"){
mask_template <- raster(nrow = 1779, ncol = 879, 
               crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ", 
               ext = extent(c(102.1417, 109.4667, 8.566667, 23.39167)))
}
if(country_abbr == "aus"){
mask_template <- raster(nrow = 4091, ncol = 4990, 
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ", 
                 ext = extent(c(112.4667, 154.05, -44.04167, -9.95)))
}

mask_template[] <- 1
mask <- mask(mask_template, mask)

if(country_abbr == "til"){
  mask <- raster(extent(90, 120, 0, 30),  
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ", 
                 res = 0.008333333)
  mask[] <- 1
}


plot(mask)

saveRDS(mask, file.path(layer_path, paste0("mask_", country_abbr, ".rds")))

#-----------------------------#
#----I. PREPARE STATIC DATA####
#-----------------------------#
mask <- readRDS(file.path(layer_path, paste0("mask_", country_abbr, ".rds")))
p
#2.a) Topography
#Data available through https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
srtm <- raster("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/srtm_aus_vnm_tile29.tif")

srtm <- projectRaster(srtm, mask)
srtm <- crop(srtm, mask)
plot(mask, add = TRUE)
srtm <- mask(srtm, mask)

names(srtm) <- "srtm"
elevation <- mask(srtm, mask)
slope <- terrain(srtm, opt = "slope")
roughness <- terrain(srtm, opt = "roughness")
terrain <- stack(elevation, slope, roughness)

#2.b) Soil
#Data availbale through https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
soil <- stack(list.files("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/soil_data/data", pattern = "*.dat", full.names = T))
crs(soil) <- crs(mask)
soil <- projectRaster(soil, mask)
soil <- crop(soil, mask)
soil <- mask(soil, mask)
names(soil) <- c( "bulk", "awco", "carb",  "nitro")

#2.c) Distance rasters and protected areas
#Shortest distance rasters calculated were distance to roads, distance to rivers, distance to built-up areas and distance to lakes.

#Data available through:
#Global road networks: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
#Global drainage systems: http://www.soest.hawaii.edu/wessel/gshhg/
#GLobal built-up areas: http://ref.data.fao.org/map?entryId=c22837d0-88fd-11da-a88f-000d939bc5d8&tab=metadata
#Global lakes: http://www.soest.hawaii.edu/wessel/gshhg/
#Protected areas: http://wcmc.io/wdpa_current_release (downloaded Feb 2018)

#Processing steps:
#1. combined shapefiles (different levels) in qgis (Levels 2-9 for rivers)
#2. rasterize in qgis (automatically subsetted from global data set)
#3. calculate proximity in qqis
#4. then load here to mask NAs and combine with rest of data

#this step is so qgis can read the attributes for some of the croppsed shapefiles
#only needs done once
# name <- "WDBII_rivers_global_L2-L9"
# country_abbr <- "vnm"
# shp <- readOGR("/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name)
# shp@data[,2] <- as.numeric(shp@data[,2])
# shp@data[,1] <- as.numeric(shp@data[,1])
# writeOGR(shp, dsn = "/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name, driver = "ESRI Shapefile", overwrite = T)

#load qgis processed files, and set NA
# lakes <- raster(paste0(raw_rasters_rivers_lakes_cl/lakes_", country_abbr, ".tif"))
# coast <- raster(paste0(raw_rasters_rivers_lakes_cl/coast_", country_abbr, ".tif"))
# rivers <- raster(paste0(raw_rasters_rivers_lakes_cl/rivers_", country_abbr, ".tif"))
# roads <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", country_abbr, ".tif"))
# builtup <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", country_abbr, ".tif"))
#names <- c("distrivers", "distlakes", "distcoastline", "PA", "distbuiltup", "distroads")

distances <- stack(list.files("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/preprocessed/", full.names = T, pattern = country_abbr)[-6])

distances <- crop(distances, mask)
distances <- mask(distances, mask)
names(distances) <- c("dibu", "dico", "dila", "diri", "diro")
rm(elevation, roughness, slope, srtm)

#PA
plot(raster(paste0("/Users/kapitzas/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/preprocessed/PA_", country_abbr, ".tif")))
path <- "/Users/kapitzas/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Global Protected Areas"
mask <- readRDS(paste0("/Users/kapitzas/OneDrive - The University of Melbourne/PhD/chapter1/birds_ccimpacts/RData/mask_", country_abbr, ".rds"))
input <- file.path(path, "WDPA_Mar2018-shapefile-polygons.shp")
output <- file.path(path, "WDPA_Mar2018-shapefile-polygons_cropped.shp")
crop_shp(input, output, extent(mask))
cropped <- readOGR(output)
subs <- cropped[which(cropped@data$IUCN_CAT%in%c("II", "Ia", "Ib")),]
output_subs <- "WDPA_Mar2018-shapefile-polygons_cropped_subs"
writeOGR(subs, path, output_subs, driver = "ESRI Shapefile", overwrite_layer = TRUE)
input <- paste0(file.path(path, output_subs), ".shp")
output <- file.path(path, "pa_raster.tif")

rasterize_shp(input, output, res = res(mask), ext = extent(mask))
ras <- raster(output)
ras <- mask(ras, mask)
writeRaster(ras, paste0("/Users/kapitzas/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/preprocessed/PA_", country_abbr, ".tif"), format = "GTiff", overwrite = TRUE)

pa <- raster(list.files("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/preprocessed/", full.names = T, pattern = country)[6])
names(pa) <- "pa"
pa <- crop(pa, mask)
pa <- mask(pa, mask)
pa[] <- (pa[] * -1) + 1

#2.d) Population density
#data available through http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
popdens <- raster("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/Pop_density/gluds00ag.bil")
popdens <- crop(popdens, mask, snap = "near")
popdens <- mask(popdens, mask)
names(popdens) <- "popd"

#2.e) Land use
#Data available through https://landcover.usgs.gov/global_climatology.php

l <- list.files(file.path("~", "OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data", "Global", "copernicus tiles", country_abbr), full.names = TRUE, recursive = FALSE)
lutiles_list <- list()
for(i in 1:length(l)){
  print(i)
  r <- raster(l[[i]])
  lutiles_list[[i]] <- projectRaster(r, mask, method = "ngb")
  names(lutiles_list[[i]]) <- "lu"
}

tiles <- tile_merge(lutiles_list)
lu <- mask(tiles, mask)

table(lutiles_list[[2]][])
lu[lu[] == 50] <- 1 #urban
lu[lu[] == 40] <- 2 #crop
lu[lu[] == 30] <- 3 #herbacious vegetation
lu[lu[] == 20] <- 4 #shurbs
lu[lu[]%in%c(121:126)] <- 5 #open forest
lu[lu[]%in%c(111:116)] <- 6 #closed forest
lu[lu[]%in%c(90, 100)] <- 7 #herbaceous wetland, moss, lichen
lu[lu[]%in%c(60, 70)] <- 8 #bare, sparese, ice
lu[lu[]%in%c(0, 80, 200)] <- NA #permanent water bodies (covered by distance to river and lake covariates)

plot(lu)


# world_lu <- raster("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/LCType.tif")
# lu <- crop(world_lu, mask)
# lu <- projectRaster(lu, mask, method = "ngb")
# lu <- mask(lu, mask)
# 
# #Reclassify, see original calsses below.
# lu[lu[]%in%c(12,14)] <- 100
# lu[lu[]%in%c(6, 7,9,10,8)] <- 200
# lu[lu[]%in%c(1,2,3,4,5)] <- 300
# lu[lu[]%in%c(13)] <- 400
# lu[lu[]%in%c(11,15,16,0)] <- 500
# 
# lu <- lu/100 - 1
# table(lu[])
names(lu) <- "landuse"
# 0	Water
# 1	Evergreen Needle leaf Forest
# 2	Evergreen Broadleaf Forest
# 3	Deciduous Needle leaf Forest
# 4	Deciduous Broadleaf Forest
# 5	Mixed Forests
# 6	Closed Shrublands
# 7	Open Shrublands
# 8	Woody Savannas
# 9	Savannas
# 10	Grasslands
# 11	Permanent Wetland
# 12	Croplands
# 13	Urban and Built-Up
# 14	Cropland/Natural Vegetation Mosaic
# 15	Snow and Ice
# 16	Barren or Sparsely Vegetated
plot(lu)
#---------------------------#
#####----II. BIOCLIM DATA####
#---------------------------#
#3. a) Download tiles, mosaic, crop and write
if(country_abbr == "aus"){
  tiles <- c("39", "310", "311", "49", "410", "411")
}else if(country_abbr%in%c("vnm", "til")){
  tiles <- "29"
}

temp_folder <- file.path(".", "temp")
dir.create(temp_folder)
wc <- get_wctiles(tile = tiles, var = "bio", path = temp_folder)
wc2 <- merge_wctiles(wc)
names(wc2) <- paste0("bio",c(1:19))
wc2 <- crop(wc2, mask)
extent(wc2) <- extent(mask)
wc2 <- mask(wc2, mask)
bioclim <- wc2
unlink(temp_folder, recursive = T)

covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)

#synch NA of present day covariates
for(i in 1:nlayers(covariates)){
  mask <- mask(mask, covariates[[i]])
  print(i)
}


for(i in 1:nlayers(covariates)){
  covariates[[i]] <- mask(covariates[[i]], mask)
  print(i)
}

summary(as.data.frame(covariates))

saveRDS(readAll(covariates), file = paste0(layer_path, "/covariates_", country_abbr, ".rds"))
saveRDS(mask, file = paste0(layer_path, "/mask_", country_abbr, ".rds"))

#3.b) Downlaoad GCM model predictions
temp_folder <- file.path(".", "temp")
dir.create(temp_folder)
regions <- c("aus", "vnm")
i <- 1
for(i in 1:length(regions)){
  assign(paste0("mask_", regions[i]), readRDS(file.path(".", "RData", paste0("mask_", regions[i], ".rds"))))
}

res <- 0.5
scens <- c("26", "85")
t <- readRDS("/Users/kapitzas/OneDrive - The University of Melbourne/PhD/chapter1/birds_ccimpacts/RData/bioq2_85_aus.rds")
plot(t)

models <- c("GF", "MP", "BC", "CC", "CN", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
download_path <- temp_folder

dir.create(download_path)
for(i in 1:length(models)){
  aus_stack <- list()
  vnm_stack <- list()
  for (j in 1:2){
    f <- tryCatch(getData("CMIP5", rcp = scens[j], year = 70, model = models[i], res = res, var = "bio", path = download_path), error = function (e) NA)
    if(class(f) != "RasterStack"){next}
    aus_stack[[j]] <- mask(crop(f, mask_aus), mask_aus)
    vnm_stack[[j]] <- mask(crop(f, mask_vnm), mask_vnm)
  }
  save(vnm_stack, file = paste0(download_path, "/vnm_", models[i], ".rda"))
  save(aus_stack, file = paste0(download_path, "/aus_", models[i], ".rda"))
  unlink(paste0(download_path, "/cmip5"), recursive = TRUE)
  rm(vnm_stack, aus_stack)
}

#ii) Extract cell-wise quartiles across GCM
country_abbr <- "aus"
gcm <- list.files("~/Dropbox/PhD - Large Files/PhD - Raw Data/GCM/", full.names = T, pattern = country_abbr)
mask <- readRDS(file.path(".", "RData", paste0("mask_", country_abbr, ".rds")))

r <- mask
inds <- which(!is.na(r[]))
quartiles <- c("q1", "q2", "q3")
scens <- c("26", "85")

k <- 1
for(k in 1:2){
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q1_", scens[k], "_", country_abbr,  ".rds"))
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q2_", scens[k], "_", country_abbr,  ".rds"))
  saveRDS(stack(), file = paste0(layer_path, "/bio", "q3_", scens[k], "_", country_abbr,  ".rds"))
  print(paste0("processing rcp", scens[k]))
  for(j in 1:19){
    print(paste0("processing cov: ", j))
    bio <- stack()
    for(i in 1:14){
      print(paste0("processing scen: ", i))
      dat <- get(load(gcm[[i]]))[[k]]
      bio <- stack(bio, dat[[j]])
    }
    
    print(paste0("getting quantiles..."))
    df1 <- na.omit(as.matrix(getValues(bio)))
    c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
    for(m in 1:3){
      bioclim <- readRDS(file = paste0(layer_path, "/bio", quartiles[m], "_", scens[k], "_", country_abbr,  ".rds"))
      r[inds] <- c[,m]
      names(r) <- paste0("bio", j)
      saveRDS(stack(bioclim, r), file = paste0(layer_path, "/bio", quartiles[m], "_", scens[k], "_", country_abbr,  ".rds"))
    }
  }
}


#------------------------------------------#
#####--III. SYNCING NA ACROSS ALL LAYERS####
#------------------------------------------#
nas <- list()
files <- list.files(layer_path, full.names = TRUE, pattern = country_abbr)
stack <- stack(files[grepl(paste0("(?=.*",country_abbr,")"), files, perl = TRUE)])
sums <- calc(stack, sum)
sums_nona <- which(is.na(sums[]))

for(i in 1:nlayers(stack)){
  r <- stack[[i]]
  r[sums_nona] <- NA
  writeRaster(r, paste0(layer_path,names(stack)[[i]]), format = "GTiff", overwrite = T)
  print(i)
}
gc()

#------------------------------------------------------#
#####--IV. MAKING BIOREGIONS RASTER (AUSTRALIA ONLY)####
#------------------------------------------------------#
#Data available through http://www.environment.gov.au/fed/catalog/search/resource/downloadData.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D
bioreg <- readOGR(file.path("~", "Dropbox", "PhD - Large Files", "PhD - Raw Data", "Australia", "IBRA7_regions", "ibra7_regions.shp"))
mask <- readRDS(file = file.path(".", "RData", "mask_aus.rds"))
bioreg_rast <- mask
bioreg_rast[] <- NA

l <- length(bioreg@polygons)
for(i in 1:l){
  print(i)
  m <- raster::extract(mask, bioreg[i,], cellnumbers = T)
  if(!is.null(m[[1]])){
  bioreg_rast[m[[1]][,1]] <- i
  }
}
bioreg_rast <- mask(bioreg_rast, mask)

saveRDS(bioreg_rast, file.path(".", "RData", "bioregions_aus.rds"))


#------------------#
#####----V. GTAP####
#------------------#
#The input data are the outputs of our CGE model
countries <- c("aus", "vnm")
rcps <- c("1oC", "4oC")
scens <- c("26", "85")
timesteps <- c(1, 13, 23, 33, 43, 53) #(#years from 2018)


#6.a) Process data for plots (2018 and 2070 sector output and landendowments)
all_bc <- data.frame(row.names = 1:58)
countries <- c("aus", "vnm")
for (l in c(1:2)){
  country_abbr <- countries[l] #"vnm"
  if (country_abbr == "aus"){
    cn <- "aus"
  }else if (country_abbr == "vnm"){
    cn <-"vn"
  }
  rcps <- c("1oC", "4oC")
  scens <- c("26", "85")
  n <- data.frame(row.names = 1:58)
  for (i in 1:length(rcps)){
    qfe <- read.table(file.path("~", "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data", "demand calcs", "prod and output", paste0("qfe_", cn, "_", rcps[i], "_all.csv")), header = F, sep = ",")
    qfe_df <- as.table(t(qfe[which(qfe$V1 == "land"),][-1]))
    coms <- t(t(qfe[2,][1,]))
    colnames(qfe_df) <- qfe$V1[seq(1,98, 8)]
    yrs <- 2070
    m <- matrix(NA, nrow(qfe_df), length(yrs))
    for(k in 1:nrow(qfe_df)){
      m[k,] <- approx(y = qfe_df[k,-13], x = as.numeric(colnames(qfe_df))[-13], xout = yrs)$y
    }
    n <- cbind(n,m)
  }
  qfe_all <- n
  rownames(qfe_all) <- coms[-1]
  
  qfe_all <- qfe_all[-which(rownames(qfe_all) == "cgds"),]
  for (i in 1:length(rcps)){
    dat <- read.table(file.path("~", "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data", "demand calcs", "prod and output", paste0("qo_",cn ,"_",rcps[i], ".csv")), sep = ",")
    coms <- trimws(as.character(dat[,2]))
    dat <- dat[c(2, na.omit(match(rownames(qfe_all), coms))), -c(1,2,3,4)]
    colnames(dat) <- dat[1,]
    dat <- dat[-1,]
    rownames(dat) <- rownames(qfe_all)
    m <- matrix(NA, nrow(dat), length(yrs))
    for(k in 1:nrow(dat)){
      m[k,] <- approx(y = dat[k,], x = as.numeric(colnames(dat)), xout = yrs)$y
    }
    m <- as.data.frame(m)
    rownames(m) <- rownames(dat)
    qfe_all <- cbind(qfe_all, m)
  }
  all <- qfe_all
  if(country_abbr == "vnm"){
    all[which(rownames(all) == "wht"),] <- 0
  }
  colnames(all) <- c(paste0(country_abbr, "_", c(paste0("qfe_",rcps), paste0("qo_", rcps))))
  all_bc <- c(all_bc, all)
}

all <- as.data.frame(all_bc, row.names = rownames(dat))
saveRDS(all, file = file.path(".", "RData", paste0("results_gtap_qo_qfe", ".rds")))

#6.b) Landendowment trajectories (for land use demand trajectories)
harvested <- read.csv("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/demand calcs/FAOSTAT_data_9-18-2018.csv", header = T)
j <- 2
i <- 2
for (j in 1:length(countries)){
  country <- countries[j]
  landuse <- readRDS(file.path(".", "RData", paste0("covariates_", country, ".rds")))
  landuse <- landuse[[which(names(landuse) == "landuse")]]
  demand_2005 <- table(landuse[])
  
  for (i in 1:length(rcps)){
    rcp <- rcps[i]
    yrs <- length(2018:2070)
    
    if (country == "aus"){
      harv <- harvested[which(harvested$Area == "Australia"),]
    }else if (country == "vnm"){
      harv <- harvested[which(harvested$Area == "Viet Nam"),]
    }
    
    #Sum the area (in ha) produced in available GTAP classes
    total_bysector <- aggregate(harv$Value, by = list(harv$GTAP.sector), FUN = function(x) sum(na.omit(x)))
    
    #relative area harvested of agricultural commodities as share of total area harvested in one of the GTAP sectors.
    
    total_bysector[,2] <- total_bysector[,2]/sum(total_bysector[,2])
    
    saveRDS(total_bysector, file.path(".", "RData", paste0("harvested2016_", country, ".rds")))
    if (country == "aus"){
      cn <- "aus"
    }else if (country == "vnm"){
      cn <-"vn"
    }
    c <- read.table(file.path("~", "OneDrive - The University of Melbourne", "PhD - Large Files", "PhD - Raw Data", "demand calcs", "prod and output", paste0("qfe_", cn, "_", rcps[i], "_all.csv")), header = F, sep = ",")
    c_df <- as.table(t(c[which(c$V1 == "land"),][-1]))
    coms <- c[2,]
    coms <- t(t(coms[1,]))
    colnames(c_df) <- c$V1[seq(1,98, 8)]
    yrs <- c(2019:2070)
    m <- matrix(NA, nrow(c_df), length(yrs))
    for(k in 1:nrow(c_df)){
      m[k,] <- approx(y = c_df[k,-13], x = as.numeric(colnames(c_df))[-13], xout = yrs)$y
    }
    
    m <- as.data.frame(cbind("2018" = 0, m))
    rownames(m) <- coms[-1]
    colnames(m)[-1] <- as.character(yrs)
    m_agr <- m[rownames(m)%in%total_bysector$Group.1,]
    m_agr <- m_agr[match(rownames(m_agr), total_bysector[,1]),]
    m_agr <- rbind(m_agr, m[which(rownames(m)%in%c("frs", "ctl")),])
    print(m_agr[,ncol(m_agr)])
    saveRDS(m_agr, file = file.path(".", "RData", paste0("gtap_landendowments_", country,"_", scens[i], ".rds")))
  }
}

#---------------------------------#
#--VI. THINNING OCCURRENCE DATA####
#---------------------------------#

#GBIF data are available via DOIs:
#Australia:  https://doi.org/10.15468/dl.khlzmu
#Tile 29: https://doi.org/10.15468/dl.yapqxq
#Viet Nam: https://doi.org/10.15468/dl.nt0ftl
#1) GLOBAL PARAMETERS####
country_abbr <- "aus" #til or aus
# 1.a) Select country and load GBIF data
occ <- as.data.frame(fread(paste0("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/GBIF/gbif_aves_", country_abbr, ".csv"), header = T, select = c("decimallongitude", "decimallatitude", "species"), na.strings=c("NA", "", " ")))
if(country_abbr == "til"){
  occ_vnm <- fread(paste0("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/GBIF/gbif_aves_vnm.csv"), header = T, select = c("decimallongitude", "decimallatitude", "species", na.strings=c("NA", "", " ")))
}

mask <- readRDS(file.path(layer_path, paste0("mask_", country_abbr, ".rds")))
#1.b) Find exact spatial duplicates and remove
species <- unique(occ$species)
occ2 <- data.frame()

#Remove dups within same species
for (i in 1:length(species)) {
  spec_inds <- which(occ$species == species[i])
  values <- extract(mask, occ[spec_inds, c(1,2)], cellnumbers = T)
  occ2 <- rbind(occ2, occ[spec_inds[which(!duplicated(values[,1]))],])
  print(i)
}

#1.c) Remove entries with less than 2 dec places (because their precision is ca. below 1km2)
lis <- strsplit(sub('0+$', '', as.character(occ2$decimallongitude)), ".", fixed = TRUE)
out <- sapply(lis, function(x) nchar(x[2]))
lis <- strsplit(sub('0+$', '', as.character(occ2$decimallatitude)), ".", fixed = TRUE)
out2 <- sapply(lis, function(x) nchar(x[2]))
unprecise <- which(out < 2 | out2 < 2)

if(length(unprecise) > 0) { occ2 <- occ2[-unprecise,] }

#1.d) check that all obs are in locations with data and remove obs where that isn't the case
bgp <- SpatialPoints(occ2[,c(1,2)])
preds_obs <- extract(mask, bgp)
occ2 <- occ2[-which(is.na(preds_obs)),]

#1.e) Filter for >= 20 occurrences
occ3 <- occ2[occ2$species%in%names(which(table(occ2$species) >= 20)),]
min(table(occ3$species)) >= 20
species_left <- unique(occ3$species)
spnums <- table(occ3$species)

dat <- occ3
names(dat) <- c("long","lat","species") #column names

#1.f) Subset the ones only occuring in vnm
if(country_abbr == "til"){
  specs <- unique(dat$species)
  specs <- specs[which(specs%in%unique(occ_vnm$species))]
  dat <- dat[dat$species%in%specs,]
}

saveRDS(dat, file = paste0(layer_path, "/occ_", country_abbr, ".rds"))


#synch NA of present day covariates
source(file.path(".", "R", "1_functions.R"))
country_abbr <- "aus" #this script needs to be run for aus and vnm
layer_path <- file.path(".", "RData")

mask <- readRDS(file.path(layer_path, paste0("mask_", country_abbr, ".rds")))
covariates <- readRDS(file.path(layer_path, paste0("covariates_", country_abbr, ".rds")))
mask <- covariates[[1]]

files_count  <- list.files(file.path(layer_path), pattern = country_abbr, full.names = TRUE)
bio <- files_count[grepl("bioq", files_count)]
covs <- files_count[grepl("covariates", files_count)]
files <- c(bio, covs)
for(j in 1:length(files)){
  r <- readRDS(files[[j]])
  for(i in 1:nlayers(r)){
    mask <- mask(mask, r[[i]])
    print(i)
  }
  print(j)
}

for(j in 1:length(files)){
  r <- readRDS(files[[j]])
  for(i in 1:nlayers(r)){
    r[[i]] <- mask(r[[i]], mask)
    print(i)
  }
  saveRDS(readAll(r), file = files[[j]])
  print(j)
}
length(which(is.na(mask[])))
length(which(is.na(r[[1]][])))
which(is.na(mask[]) != is.na(r[[1]][]))

