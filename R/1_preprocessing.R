## Processing input data for SDMs
##
## Outputs:
##  1. masks (til, vn, aus) - raster layer [3]
##  2. covariates (til, vn, aus) - raster stack [3]
##  3. bioq (0.25, 0.5 and 0.75 quartiles for vn, aus for each rcp) 
##      - raster stack [3 x 2 x 3 = 18]
##  4. bioregions (aus) - raster layer [1]
##  5. harvested2016 (vn, aus) - data.frame (prop of land harvested for 
##      commodity x in 2016) [8 commodities x 2]
##  6. gtap_landendowments (vn, aus for rcp 26 & 85) - data.frame 
##      (prop of land harvested 'endowmnet for x commodities for each year 
##      from 2019 - 2070) [10 sectors/commodities x 53 years]
##  7. occ_ (vn, aus) - bird data from gbif - data.frame


## Set up work environment ####
# setwd("./regtrade/")
rm(list = ls())
# devtools::install_github('kapitzas/WorldClimTiles')
# devtools::install_github('skiptoniam/sense')
# devtools::install_github('smwindecker/gdaltools')
x <- c('data.table','rgdal','rgeos','matrixStats','rgdal',"sp",'raster',
       'WorldClimTiles','sense', 'gdaltools' , 'readxl')
lapply(x, require, character.only = TRUE)
source(file.path(".", "R", "1_functions.R"))
gsdms_data <- "/Volumes/discovery_data/gsdms_data" # change as per server: boab - "./data"
regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' # change as per server: boab = "./data"
layer_path <- file.path(regSSP_birds_data, "RData") # change as per server: boab - "./RDdata"
if(!dir.exists(layer_path)){dir.create(layer_path)}

## 1. Prepare masks ####
regions <- c("vn", "aus", "til")

for (region in regions){
  
  if (region == "til") {
    reg_mask <- raster(extent(90, 120, 0, 30),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
                       res = 0.008333333)
    reg_mask[] <- 1
    saveRDS(reg_mask, file.path(layer_path, paste0("mask_", region, ".rds")))
  } else {
    reg_mask <- getData("GADM", country = region, level = 0, path = layer_path)
    reg_mask <- gSimplify(reg_mask, tol = 0.00833)
    
    if (region == "vn") {
      mask_template <- raster(
        nrow = 1779,
        ncol = 879,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
        ext = extent(c(102.1417, 109.4667, 8.566667, 23.39167))
      )
      mask_template[] <- 1
    }
    
    else if (region == "aus") {
      mask_template <- raster(
        nrow = 4091,
        ncol = 4990,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
        ext = extent(c(112.4667, 154.05,-44.04167,-9.95))
      )
      mask_template[] <- 1
    }
    reg_mask <- mask(mask_template, reg_mask)
    saveRDS(reg_mask, file.path(layer_path, paste0("mask_", region, ".rds")))
    plot(reg_mask)
  }
}
rm(mask_template)


## ***** Specify region for analysis (up till Bioclim - future) ***** ####
##  Rerun for c("vn", "aus", "til")
region <- 'vn'


## 2. Prepare covariate data ####
reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))

## 2a. Topography ####
## Source: https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
## Select polygon including vn and aus
srtm <- raster(file.path(regSSP_birds_data, "sdat_10008_1_20191110_184812575.asc"))
crs(srtm) <- crs(reg_mask)
srtm <- projectRaster(srtm, reg_mask)
srtm <- crop(srtm, reg_mask)
srtm <- mask(srtm, reg_mask)
plot(srtm)

names(srtm) <- "srtm"
elevation <- mask(srtm, reg_mask)
slope <- terrain(srtm, opt = "slope")
roughness <- terrain(srtm, opt = "roughness")
terrain <- stack(elevation, slope, roughness)
rm(elevation, roughness, slope, srtm)

## 2b. Soil ####
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
## See gsdms for data download instructions
soil <- stack(list.files(file.path(gsdms_data, "orders"), pattern = "*.dat", full.names = T, recursive = T))
crs(soil) <- crs(reg_mask)
soil <- projectRaster(soil, reg_mask)
soil <- crop(soil, reg_mask)
soil <- mask(soil, reg_mask)
names(soil) <- c( "bulk", "awco", "carb",  "nitro")


## 2c. Distance rasters ####
## Shortest distance rasters calculated were distance to roads, distance to rivers, distance to built-up areas and distance to lakes.

## Sources:
## Global road networks: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
## Global drainage systems: http://www.soest.hawaii.edu/wessel/gshhg/
## GLobal built-up areas: http://ref.data.fao.org/map?entryId=c22837d0-88fd-11da-a88f-000d939bc5d8&tab=metadata
## Global lakes: http://www.soest.hawaii.edu/wessel/gshhg/
## Protected areas: http://wcmc.io/wdpa_current_release (downloaded Feb 2018)

  ## Processing steps:
  ## 1. combined shapefiles (different levels) in qgis (Levels 2-9 for rivers)
  ## 2. rasterize in qgis (automatically subsetted from global data set)
  ## 3. calculate proximity in qqis
  ## 4. then load here to mask NAs and combine with rest of data
  
  ## this step is so qgis can read the attributes for some of the croppsed shapefiles
  ## only needs done once
  # name <- "WDBII_rivers_global_L2-L9"
  # region <- "vn"
  # shp <- readOGR("/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name)
  # shp@data[,2] <- as.numeric(shp@data[,2])
  # shp@data[,1] <- as.numeric(shp@data[,1])
  # writeOGR(shp, dsn = "/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name, driver = "ESRI Shapefile", overwrite = T)
  ## load qgis processed files, and set NA
  # lakes <- raster(paste0(raw_rasters_rivers_lakes_cl/lakes_", region, ".tif"))
  # coast <- raster(paste0(raw_rasters_rivers_lakes_cl/coast_", region, ".tif"))
  # rivers <- raster(paste0(raw_rasters_rivers_lakes_cl/rivers_", region, ".tif"))
  # roads <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", region, ".tif"))
  # builtup <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", region, ".tif"))
  # names <- c("distrivers", "distlakes", "distcoastline", "PA", "distbuiltup", "distroads")

distances <- stack(list.files(file.path(regSSP_birds_data, "qgis_files"), full.names = T, pattern = region)[-6])
distances <- crop(distances, reg_mask)
distances <- mask(distances, reg_mask)
names(distances) <- c("dibu", "dico", "dila", "diri", "diro")


## 2d. Protected areas #### - to be fixed
# file_in <- list.files(file.path(gsdms_data, "protectedareas"), pattern= "WDPA_Nov2018-shapefile-polygons.shp", full.names = T, recursive = T)
# file_out <- file.path(regSSP_birds_data, "qgis_files", basename(file_in))
# file_out <- gsub(".shp", "_cropped.shp", file_out)
# crop_shp(file_in, file_out, extent(reg_mask))
# cropped <- readOGR(file_out)
# subs <- cropped[which(cropped@data$IUCN_CAT%in%c("II", "Ia", "Ib")),]
# output_subs <- "WDPA_Mar2018-shapefile-polygons_cropped_subs"
# writeOGR(subs, file.path(regSSP_birds_data, "qgis_files"), output_subs, driver = "ESRI Shapefile", overwrite_layer = FALSE)
# 
# file_in <- paste0(file.path(regSSP_birds_data, "qgis_files", output_subs), ".shp")
# file_out <- file.path(regSSP_birds_data, "qgis_files", paste0("pa_raster_", region, ".tif"))
# gdaltools::rasterize_shp(file_in, file_out, res = res(reg_mask), ext = extent(reg_mask))
# 
# pa <- raster(file_out)
# pa <- crop(pa, reg_mask)
# pa <- mask(pa, reg_mask)
# pa[] <- (pa[] * -1) + 1
# names(pa) <- "pa"
# # writeRaster(pa, paste0(layer_path, "/pa_", region, ".tif"), format = "GTiff", overwrite = TRUE)
# rm(cropped, subs, output_subs)
pa <- list.files(file.path(regSSP_birds_data, "qgis_files"), full.names = T, pattern = region)
pa <- pa[grepl("PA_", pa)]
pa <- raster(pa)
names(pa) <- "pa"
pa <- crop(pa, reg_mask)
pa <- mask(pa, reg_mask)
pa[] <- (pa[] * -1) + 1



## 2e. Population density ####
## Source: http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
popdens <- raster(file.path(gsdms_data, "popdens", "gluds00ag.bil"))
popdens <- crop(popdens, reg_mask, snap = "near")
popdens <- mask(popdens, reg_mask)
names(popdens) <- "popd"


## 2f. Land use ####
## Source: https://archive.usgs.gov/archive/sites/landcover.usgs.gov/global_climatology.html
## direct link: https://archive.usgs.gov/archive/sites/landcover.usgs.gov/documents/GlobalLandCover_tif.zip
l <- list.files(file.path(regSSP_birds_data, "copernicus_tiles", region), full.names = TRUE, recursive = FALSE)
lutiles_list <- list()
for(i in 1:length(l)){
  print(i)
  r <- raster(l[[i]])
  lutiles_list[[i]] <- projectRaster(r, reg_mask, method = "ngb")
  names(lutiles_list[[i]]) <- "lu"
}

tiles <- tile_merge(lutiles_list)
lu <- mask(tiles, reg_mask)

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
names(lu) <- "landuse"
plot(lu)
rm(r, lutiles_list, tiles)

# ## Source: https://earthexplorer.usgs.gov/
# ## See gsdms script for dowload instructions 
# world_lu <- raster("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/LCType.tif")
# lu <- crop(world_lu, reg_mask)
# lu <- projectRaster(lu, reg_mask, method = "ngb")
# lu <- mask(lu, reg_mask)
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
# names(lu) <- "landuse"
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


## 2g. Bioclim - current ####
## Download tiles, mosaic, crop and write
file_in <- list.files(file.path(gsdms_data, 'bio_30s'), full.names = T)
bioclim <- list()

for (f in 1:length(file_in)){
  bioclim[[f]] <- crop(raster(file_in[f]), reg_mask)
  bioclim[[f]] <- mask(bioclim[[f]], reg_mask)
}
bioclim <- stack(bioclim)
  

## Sync NAs - I ####
## Find min non-NA set values across mask and covariates and sync NAs 
covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)
for(i in 1:nlayers(covariates)){
  reg_mask <- mask(reg_mask, covariates[[i]]) # to find minimum non-NA set values for mask layer
  print(i)
}

for(i in 1:nlayers(covariates)){
  covariates[[i]] <- mask(covariates[[i]], reg_mask) # mask using updated mask layer
  print(i)
}

summary(reg_mask) # summary(reg_mask[])
summary(as.data.frame(covariates))
saveRDS(readAll(covariates), file = paste0(layer_path, "/covariates_", region, ".rds"))
saveRDS(reg_mask, file = paste0(layer_path, "/mask_", region, ".rds"))

rm(terrain, soil, distances, pa, popdens, bioclim, lu)


## ***** Rerun scipt up till here for all regions ***** ####


## 2h. Bioclim - future ####

## Global variables
regions <- c("vn", "aus")
rcps <- c("45", "60", "85")
models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")


  # ## Downlaoad GCM model predictions
  # urls_bioclim <- c()
  # for (model_i in models){
  #   for (rcp_j in rcps){
  #     fn = paste0(tolower(model_i), rcp_j, 'bi70', '.zip')
  #     thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/cmip5/30s/', fn)
  #     urls_bioclim = c(urls_bioclim, thisurl)
  #   }
  # }
  # library(bitops)
  # library(RCurl)
  # test = lapply(urls_bioclim, url.exists)
  # all(unlist(test))
  # 
  # zipdst = file.path(gsdms_data, 'gcm_30s.zip')
  # rasterdst = file.path(gsdms_data, 'gcm_30s/')
  # if(!dir.exists(rasterdst)) {
  #   dir.create(rasterdst)
  # } # end if !dir.exists
  # for (i in 1:length(urls_bioclim)){
  #   download_from_url(urls = urls_bioclim[i], zipdst, rasterdst)
  # }
  # 
  # ## Mask data for study regions & create stacks by region and gcms
  # files_gcm <- list.files(file.path(gsdms_data, 'gcm_30s'), full.names = T, recursive = T)
  # gcm_reg_path <- file.path(regSSP_birds_data, 'gcm_reg')
  # dir.create(gcm_reg_path)
  # 
  # for(region in regions){
  #   for(model_i in models){
  #     reg_stack <- list()
  #     file_mod <- files_gcm[grepl(model_i, files_gcm)]
  #     for(j in 1:length(rcps)){
  #       file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
  #       temp_stack <- list()
  #       for(f in 1:length(file_mod_rcp)){
  #         reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
  #         temp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), reg_mask), reg_mask)
  #       }
  #       reg_stack[[j]] <- readAll(brick(temp_stack))
  #     }
  #     saveRDS(reg_stack, file = paste0(gcm_reg_path, "/", region, "_", model_i, ".rds"))
  #   }
  # }
  # rm(temp_stack, reg_stack)
  # 
  # ## Extract cell-wise quartiles across GCM
  # quartiles <- c("q1", "q2", "q3")
  # 
  # for(region in regions){
  #   gcm <- list.files(gcm_reg_path, full.names = T, pattern = region)
  #   reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
  #   r <- reg_mask
  #   inds <- which(!is.na(r[]))
  #   
  #   j<-1
  #   for(j in 1:length(rcps)){
  #     saveRDS(stack(), file = paste0(layer_path, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
  #     saveRDS(stack(), file = paste0(layer_path, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
  #     saveRDS(stack(), file = paste0(layer_path, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
  #     print(paste0("processing rcp", rcps[j]))
  #     for(k in 1:19){
  #       print(paste0("processing bioclim var: ", k))
  #       bio <- stack()
  #       for(i in 1:length(models)){
  #         print(paste0("processing model: ", i))
  #         dat <- readRDS(gcm[[i]])[[j]]
  #         bio <- stack(bio, dat[[k]])
  #       }
  #       
  #       print(paste0("getting quartiles..."))
  #       df1 <- na.omit(as.matrix(getValues(bio)))
  #       c <-rowQuartiles(df1, probs = c(0.25, 0.5, 0.75))
  #       for(m in 1:3){
  #         bioclim <- readRDS(file = paste0(layer_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
  #         r[inds] <- c[,m]
  #         names(r) <- paste0("bio", k)
  #         saveRDS(readAll(stack(bioclim, r)), file = paste0(layer_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
  #       }
  #     }
  #   }
  # }
  # # unlink(gcm_reg_path, recursive=T)


## Sync NAs - II ### 
## For covariartes + bioq layers
## Find min non-NA set values across mask, covariates and bioq layers and sync NAs 
regions <- c("vn", "aus") ## or c("vn", "aus", "til")
for (region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
  covariates <- readRDS(file.path(layer_path, paste0("covariates_", region, ".rds")))
  
  files  <- list.files(file.path(layer_path), pattern = region, full.names = TRUE)
  bioq <- files[grepl("bioq", files)]
  covs <- files[grepl("covariates", files)]
  files <- c(bioq, covs)
  for(j in 1:length(files)){
    print(paste0("processing file = ", j, " of ", length(files) ,":", files[j]))
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      print(paste0("processing layer = ", i, " of ", nlayers(r)))
      reg_mask <- mask(reg_mask, r[[i]])
    }
  }
  saveRDS(reg_mask, file = paste0(layer_path, "/mask_", region, ".rds"))
  
  for(j in 1:length(files)){
    print(paste0("processing file = ", j, " of ", length(files) ,":", files[j]))
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      print(paste0("processing layer = ", i, " of ", nlayers(r)))
      r[[i]] <- mask(r[[i]], reg_mask)
    }
    saveRDS(readAll(r), file = files[[j]])
  }
  print(paste0("summary for ", region, "... "))
  print(paste0("NAs in mask = ", length(which(is.na(reg_mask[])))))
  print(paste0("NAs in stack = ", length(which(is.na(r[[1]][])))))
  ## only check againt first layer of stack r[[1]] (all layers inn stack are the same)
  print(paste0("Which NAs in mask != NAs in layers: ", length(which(is.na(reg_mask[]) != is.na(r[[1]][])))))
  print("=========================")
}

summary(reg_mask) # summary(reg_mask[])
summary(as.data.frame(covariates))

rm(files, bioq, covs, i, j, r, regions, reg_mask, covariates)

## Reduce raster size ####
##  Remove NA values from covariate stack 
##  Mask retained with NA values
regions <- c("vn", "aus", "til")
layer_path2 <- file.path(regSSP_birds_data, "nonatables") # change as per server: boab - "./nonatables"
if(!dir.exists(layer_path2)){dir.create(layer_path2)}

for(region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
  reg_layers <- list.files(layer_path, full.names = TRUE, pattern = region)
  ## files[grepl(paste0("(?=.*",region,")"), files, perl = TRUE)]
  reg_layers <- reg_layers[!grepl("mask", reg_layers)]
  
  for(j in 1:length(reg_layers)){
    print(paste0("processing file = ", j, " of ", length(reg_layers) ,":", reg_layers[j]))
    reg_stack <- readRDS(reg_layers[j])
    ind_nona <- which(!is.na(reg_mask[]))
    
    print(paste0("# layers in ", reg_layers[j], " = ", nlayers(reg_stack)))
    new_dat <- data.table()
    for(i in 1:nlayers(reg_stack)){
      print(paste0("processing layer = ", i, " of ", nlayers(reg_stack)))
      r <- reg_stack[[i]]
      r <- getValues(r)
      r <- r[ind_nona]
      new_dat[,paste0("new_dat", "$", names(reg_stack[[i]])) := r]    
    }
    saveRDS(new_dat, file = paste0(layer_path2, basename(reg_layers[j]), "_nona.", "rds"))
    ## txt or csv files using fwrite much larger!
  }
}
rm(layer_path2, reg_mask, reg_layers, reg_stack, ind_nona, nnew_dat, r)
gc()



## Save raster stacks as data.table ####



## 3. Bioregions layer (for Australia only) ####
## Source: http://www.environment.gov.au/fed/catalog/search/resource/downloadData.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D
bioreg <- readOGR(file.path(regSSP_birds_data, "IBRA7_regions", "ibra7_regions.shp"))
reg_mask <- readRDS(file = file.path(layer_path, "mask_aus.rds"))
bioreg_rast <- reg_mask
bioreg_rast[] <- NA

l <- length(bioreg@polygons)
for(i in 1:l){
  print(i)
  m <- raster::extract(reg_mask, bioreg[i,], cellnumbers = T)
  if(!is.null(m[[1]])){
    bioreg_rast[m[[1]][,1]] <- i
  }
}
bioreg_rast <- mask(bioreg_rast, reg_mask)
saveRDS(bioreg_rast, file.path(layer_path, "bioregions_aus.rds"))


## 4. GTAP data ####
regions <- c("aus", "vn")
ssps <- paste0("ssp", 1:3)

## Area harvested for crops in 2016 - FAO data
## Data source:http://www.fao.org/faostat/en/#data/QC
## 'Item' in data is coded for GTAP sectors/commodities befor importing in R
harvested <- read.csv(file.path(regSSP_birds_data,"fao_data/FAOSTAT_data_9-18-2018.csv"), header = T) 

for (region in regions){
  if (region == "aus"){
    harv <- harvested[which(harvested$Area == "Australia"),]
  }else if (region == "vn"){
    harv <- harvested[which(harvested$Area == "Viet Nam"),]
  }
  
  ## Estimate sum of area harvested (in ha) for each GTAP sector
  total_bysector <- aggregate(harv$Value, by = list(harv$GTAP.sector), FUN = function(x) sum(na.omit(x)))
  
  ## Estimate relative area harvested (per GTAP sector) relative to total area harvested for all GTAP sectors considered (7 sectors)
  total_bysector[,2] <- total_bysector[,2]/sum(total_bysector[,2])
  colnames(total_bysector) <- c("gtap_sector", "prop_harvested")
  saveRDS(total_bysector, file.path(layer_path, paste0("harvested2016_", region, ".rds")))
  
  ## Specify additional gtap sectors (i.e. forestry, cattle)
  gtap_sectors <- c(as.character(total_bysector$gtap_sector), "frs", "ctl")
  
  ## Load GTAP land endowmenet by sector data for 2020 - 2071
  ## Source: Tom Kompas (Uni Melb), Van Ha Pham (ANU)
  for (ssp in ssps){
    temp <- as.data.table(read_xls(file.path(regSSP_birds_data, "gtap_data", "ssps",
                                             paste0(region, "_", ssp, ".xls")), sheet = "land"))
    temp <- temp[X__1 %in% gtap_sectors]
    temp[, c(2:6, ncol(temp)) := NULL]
    temp[, c(20:29) := NULL]
    colnames(temp)[1] <- "gtap_sector"
  
    ## Estimate land endowment values for each year from 2020 - 2070 using linear interpolation
    new_yrs <- c(2020:2070)
    new_temp <- matrix(NA, nrow(temp), length(new_yrs))
    for(k in 1:nrow(temp)){
      new_temp[k,] <- approx(y = temp[k,-1], x = as.numeric(colnames(temp)[-1]), xout = new_yrs)$y
    }
    rownames(new_temp) <- temp$gtap_sector
    colnames(new_temp) <- as.character(new_yrs)
    
    ## Add base year column for 2019
    new_temp <- cbind("2019" = 0, new_temp)

    ## Save output data
    saveRDS(new_temp, file = file.path(layer_path, paste0("gtap_landendowments_", region,"_", ssp, ".rds")))
    
    ## print land endowments for 2070
    print(paste0("land endowments for 2070 for ", region, " and ", ssp, ":"))
    print(new_temp[,ncol(new_temp)])
    print("----------------------") 
  }
}


## 7. Biodiversity data ####
## Australia:  https://doi.org/10.15468/dl.khlzmu
## Tile 29: https://doi.org/10.15468/dl.yapqxq
## Viet Nam: https://doi.org/10.15468/dl.nt0ftl
source('./R/filter.gbifcsv.R')
gbif_backbone <- fread("/Volumes/discovery_data/discovery_paper_1/data/raw/gbif/gbif_backbone_taxonomy.tsv")

for(region in regions){
  if(region == "til"){region <- 'vn'}
  
  gbif_all <- fread(paste0(file.path(regSSP_birds_data, "gbif/gbif_aves_"), region, ".csv"), header = T, na.strings=c("NA", "", " "))

filter.gbifcsv(gbif.downloaded.data = gbif_all,
               gbif.nub.taxonomy = gbif_backbone,
               subset.gbifnubtaxonomy.byclass = "Aves",
               output_folder = "./data",
               output_name = paste0("occ_", region, ".rds"),
               domain.mask = reg_mask,
               select_fields = c("decimallongitude", "decimallatitude", "species"))
}
rm(gbif_backbone)






## EXTRAS ###

## Bioclim current - Download and mask
# if(region == "aus"){
#   tiles <- c("39", "310", "311", "49", "410", "411")
# }else if(region%in%c("vn", "til")){
#   tiles <- "29"
# }
# 
# temp_folder <- file.path(".", "temp")
# dir.create(temp_folder)
# wc <- get_wctiles(tile = tiles, var = "bio", path = temp_folder)
# wc2 <- merge_wctiles(wc)
# names(wc2) <- paste0("bio",c(1:19))
# wc2 <- crop(wc2, reg_mask)
# extent(wc2) <- extent(reg_mask)
# wc2 <- mask(wc2, reg_mask)
# bioclim <- wc2
# unlink(temp_folder, recursive = T)


## Bioclim future - Mask data for study regions & stack
# dir.create("./temp")
# for(region in regions){
#   assign(paste0("mask_", region), readRDS(file.path(layer_path, paste0("mask_", region, ".rds"))))
#   
#   file_in <- list.files(file.path(gsdms_data, 'gcm_30s'), full.names = T, recursive = T)
#   file_out <- paste0("./temp/", paste0(tools::file_path_sans_ext(basename(file_in)), "_", region, ".", tools::file_ext(file_in)))
#   
#   e <- extent(get(paste0("mask_", region)))
#   reso <- res(get(paste0("mask_", region)))
#   mapply(gdal_crop, inpath = file_in, outpath = file_out, MoreArgs = list(extent=e,res=reso))  # does nto crop vietnam land mass, only crops by bounding box
#   # mapply(gdal_mask, inpath = file_out, outpath = file_out, mask = get(paste0("mask_", region))) # not working...
# }
# rm(file_in, file_out, e, reso)
# 
# ## Create raster stacks by region and gcms
# dir.create("./temp2")
# file_out <- list.files("./temp", full.names = T)
# models <-tolower(models)
# 
# for(region in regions){
#   for (model_i in models){
#     reg_stack <- list()
#     file_reg <- file_out[grepl(region, file_out)]
#     file_reg <- file_reg[grepl(model_i, file_reg)]
#     for(j in 1:length(rcps)){
#       file_reg_rcp <- file_reg[grepl(rcps[j], file_reg)]
#       reg_stack[[j]] <- stack(file_reg_rcp)
#     }
#     save(reg_stack, file = paste0("./temp2/", region, "_", model_i, ".rds"))
#   }
# }


# ## Biodiversity data
# ## Select country and load GBIF data
# gbif_dat <- as.data.frame(fread(paste0(file.path(regSSP_birds_data, "gbif/gbif_aves_"), region, ".csv"), header = T, select = c("decimallongitude", "decimallatitude", "species"), na.strings=c("NA", "", " ")))
# 
# reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
# ## Find exact spatial duplicates and remove
# species <- unique(gbif_dat$species)
# occ2 <- data.frame()
# 
# ## Remove dups within same species
# for (i in 1:length(species)) {
#   spec_inds <- which(occ$species == species[i])
#   values <- extract(reg_mask, occ[spec_inds, c(1,2)], cellnumbers = T)
#   occ2 <- rbind(occ2, occ[spec_inds[which(!duplicated(values[,1]))],])
#   print(i)
# }
# 
# ## Remove entries with less than 2 dec places (because their precision is ca. below 1km2)
# lis <- strsplit(sub('0+$', '', as.character(occ2$decimallongitude)), ".", fixed = TRUE)
# out <- sapply(lis, function(x) nchar(x[2]))
# lis <- strsplit(sub('0+$', '', as.character(occ2$decimallatitude)), ".", fixed = TRUE)
# out2 <- sapply(lis, function(x) nchar(x[2]))
# unprecise <- which(out < 2 | out2 < 2)
# 
# if(length(unprecise) > 0) { occ2 <- occ2[-unprecise,] }
# 
# ## check that all obs are in locations with data and remove obs where that isn't the case
# bgp <- SpatialPoints(occ2[,c(1,2)])
# preds_obs <- extract(reg_mask, bgp)
# occ2 <- occ2[-which(is.na(preds_obs)),]
# 
# ## Filter for >= 20 occurrences
# occ3 <- occ2[occ2$species%in%names(which(table(occ2$species) >= 20)),]
# min(table(occ3$species)) >= 20
# species_left <- unique(occ3$species)
# spnums <- table(occ3$species)
# 
# dat <- occ3
# names(dat) <- c("long","lat","species") #column names
# 
# ## Subset the ones only occuring in vn
# if(region == "til"){
#   specs <- unique(dat$species)
#   specs <- specs[which(specs%in%unique(occ$species))]
#   dat <- dat[dat$species%in%specs,]
# }
# 
# saveRDS(dat, file = paste0(layer_path, "/occ_", region, ".rds"))
