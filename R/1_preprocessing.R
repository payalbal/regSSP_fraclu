## Script for processing input data for SDMs
##
## Outputs:
##  1. masks (til, vnm, aus) - raster layer
##  2. covariates (til, vnm, aus) - raster stack 
##  3. bioq (vnm, aus for rcp 45, 6 & 85) - raster stack 
##      (19 future bioclim varaibles )
##  4. bioregions (aus) - raster layer
##  5. harvested2016 (vnm, aus) - data.frame (prop of land harvested for 
##      commodity x in 2016) [8 commodities x 2]
##  6. gtap_landendowments (vnm, aus for rcp 26 & 85) - data.frame 
##      (prop of land harvested 'endowmnet for x commodities for each year 
##      from 2019 - 2070) [10 commodities x 53 years]
##  7. occ_ (vnm, aus) - bird data from gbif 


## Set up work environment ####
rm(list = ls())
# devtools::install_github('kapitzas/WorldClimTiles')
# devtools::install_github('skiptoniam/sense')
# devtools::install_github('smwindecker/gdaltools')
x <- c('data.table','rgdal','rgeos','matrixStats','rgdal',"sp",'raster',
       'WorldClimTiles','sense', 'gdaltools')
lapply(x, require, character.only = TRUE)
source(file.path(".", "R", "1_functions.R"))
layer_path <- file.path(".", "RData")
gsdms_data <- "/Volumes/discovery_data/gsdms_data"
regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data'


## Specify global variable ####
regions <- c("aus", "vnm")
rcps <- c("45", "60", "85")


## 1. Prepare masks ####
for (region in c("til", "vnm", "aus")){
  reg_mask <- getData("GADM", country = region, level = 0, path = layer_path)
  reg_mask <- gSimplify(reg_mask, tol = 0.00833)
  
  if (region == "vnm") {
    mask_template <- raster(
      nrow = 1779,
      ncol = 879,
      crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
      ext = extent(c(102.1417, 109.4667, 8.566667, 23.39167))
    )
    mask_template[] <- 1
    reg_mask <- mask(mask_template, reg_mask)
  }
  
  if (region == "aus") {
    mask_template <- raster(
      nrow = 4091,
      ncol = 4990,
      crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
      ext = extent(c(112.4667, 154.05,-44.04167,-9.95))
    )
    mask_template[] <- 1
    reg_mask <- mask(mask_template, reg_mask)
  }
  
  if (region == "til") {
    reg_mask <- raster(extent(90, 120, 0, 30),
                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
                   res = 0.008333333)
    reg_mask[] <- 1
  }
  
  saveRDS(reg_mask, file.path(layer_path, paste0("mask_", region, ".rds")))
  plot(reg_mask)
}


## ***** Specify region for analysis (up till Bioclim - future) ***** ####
##  Rerun for c("vnm", "aus", "til")
region <- 'vnm'


## 2. Prepare covariate data ####
reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))

## 2a. Topography ####
## Source: https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
srtm <- raster(file.path(regSSP_birds_data, "srtm_aus_vnm_tile29.tif"))

srtm <- projectRaster(srtm, reg_mask)
srtm <- crop(srtm, reg_mask)
plot(reg_mask, add = TRUE)
srtm <- mask(srtm, reg_mask)

names(srtm) <- "srtm"
elevation <- mask(srtm, reg_mask)
slope <- terrain(srtm, opt = "slope")
roughness <- terrain(srtm, opt = "roughness")
terrain <- stack(elevation, slope, roughness)


## 2b. Soil ####
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
## See gsdms for data download instructions
soil <- stack(list.files(file.path(gsdms_data, "orders", pattern = "*.dat", full.names = T, recursive = T)))
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
# region <- "vnm"
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
rm(elevation, roughness, slope, srtm)


## 2d. Protected areas ####
plot(raster(paste0(regSSP_birds_data, "/qgis_files/PA_", region, ".tif")))
reg_mask <- readRDS(file.path(".", "RData", paste0("mask_", region, ".rds")))
file_in <- list.files(file.path(gsdms_data, "protectedareas"), pattern= "WDPA_Nov2018-shapefile-polygons.shp", full.names = T, recursive = T)
file_out <- gsub(".shp", "_cropped.shp", file_in)
crop_shp(file_in, file_out, extent(reg_mask))
cropped <- readOGR(file_out)
subs <- cropped[which(cropped@data$IUCN_CAT%in%c("II", "Ia", "Ib")),]
output_subs <- "WDPA_Mar2018-shapefile-polygons_cropped_subs"
writeOGR(subs, path, output_subs, driver = "ESRI Shapefile", overwrite_layer = TRUE)
file_in <- paste0(file.path(path, output_subs), ".shp")
file_out <- file.path(path, "pa_raster.tif")

gdaltools::rasterize_shp(file_in, file_out, res = res(reg_mask), ext = extent(reg_mask))
ras <- raster(file_out)
ras <- mask(ras, reg_mask)
writeRaster(ras, paste0(regSSP_birds_data, "/qgis_files/PA_", region, ".tif"), format = "GTiff", overwrite = TRUE)

pa <- raster(list.files(file.path(regSSP_birds_data, "qgis_files"), full.names = T, pattern = region)[6])
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
l <- list.files(file.path(gsdms_data, "copernicus_tiles", region), full.names = TRUE, recursive = FALSE)
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

plot(lu)

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


## 2g. Bioclim - current ####
## Download tiles, mosaic, crop and write
file_in <- list.files(file.path(gsdms_data, 'bio_30s'), full.names = T)
bioclim <- list()

for (f in 1:length(file_in)){
  bioclim[[f]] <- crop(raster(file_in[f]), reg_mask)
  bioclim[[f]] <- mask(bioclim[[f]], reg_mask)
}
bioclim <- stack(bioclim)
  
## Sync NAs
covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)

for(i in 1:nlayers(covariates)){
  reg_mask <- mask(reg_mask, covariates[[i]]) # to find minimum set non NA values for mask
  print(i)
}

for(i in 1:nlayers(covariates)){
  covariates[[i]] <- mask(covariates[[i]], reg_mask) # mask using minimum set mask
  print(i)
}

summary(as.data.frame(covariates))
saveRDS(readAll(covariates), file = paste0(layer_path, "/covariates_", region, ".rds"))
saveRDS(reg_mask, file = paste0(layer_path, "/mask_", region, ".rds"))



## ***** Rerun scipt up till here for all regions ***** ####


## 2h. Bioclim - future ####

## Global variables
rcps <- c("45", "60", "85")
models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
regions <- c("vnm", "aus")

## Downlaoad GCM model predictions
urls_bioclim <- c()
for (model_i in models){
  for (rcp_j in rcps){
    fn = paste0(tolower(model_i), rcp_j, 'bi70', '.zip')
    thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/cmip5/30s/', fn)
    urls_bioclim = c(urls_bioclim, thisurl)
  }
}
library(bitops)
library(RCurl)
test = lapply(urls_bioclim, url.exists)
all(unlist(test))

zipdst = file.path(gsdms_data, 'gcm_30s.zip')
rasterdst = file.path(gsdms_data, 'gcm_30s/')
if(!dir.exists(rasterdst)) {
  dir.create(rasterdst)
} # end if !dir.exists
for (i in 1:length(urls_bioclim)){
  download_from_url(urls = urls_bioclim[i], zipdst, rasterdst)
}

## Mask data for study regions & create stacks by region and gcms
files_gcm <- list.files(file.path('./data', 'gcm_30s'), full.names = T, recursive = T)
gcm_reg_path <- "./data/gcm_reg/"
dir.create(gcm_reg_path)

for(region in regions){
  for(model_i in models){
    reg_stack <- list()
    file_mod <- files_gcm[grepl(model_i, files_gcm)]
    for(j in 1:length(rcps)){
      file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
      temp_stack <- list()
      for(f in 1:length(file_mod_rcp)){
        reg_mask <- readRDS(file.path(".", "RData", paste0("mask_", region, ".rds")))
        temp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), reg_mask), reg_mask)
      }
      reg_stack[[j]] <- readAll(brick(temp_stack))
    }
    saveRDS(reg_stack, file = paste0(gcm_reg_path, region, "_", model_i, ".rds"))
  }
}
rm(temp_stack, reg_stack)


## Extract cell-wise quartiles across GCM
layer_path <- "./RData"
quartiles <- c("q1", "q2", "q3")

for(region in regions){
  gcm <- list.files('/Users/payalb/Dropbox/discovery_trade/analyses/regSSP_birds/temp/gcm_reg', full.names = T, pattern = region)
  reg_mask <- readRDS(file.path(".", "RData", paste0("mask_", region, ".rds")))
  r <- reg_mask
  inds <- which(!is.na(r[]))
  
  for(j in 1:length(rcps)){
    saveRDS(stack(), file = paste0(layer_path, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
    saveRDS(stack(), file = paste0(layer_path, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
    saveRDS(stack(), file = paste0(layer_path, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
    print(paste0("processing rcp", rcps[k]))
    for(k in 1:19){
      print(paste0("processing bioclim var: ", k))
      bio <- stack()
      for(i in 1:length(models)){
        print(paste0("processing model: ", i))
        dat <- readRDS(gcm[[i]])[[j]]
        bio <- stack(bio, dat[[k]])
      }
      
      print(paste0("getting quantiles..."))
      df1 <- na.omit(as.matrix(getValues(bio)))
      c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
      for(m in 1:3){
        bioclim <- readRDS(file = paste0(layer_path, "/bio", quartiles[m], "_", 
                                         rcps[j], "_", region,  ".rds"))
        r[inds] <- c[,m]
        names(r) <- paste0("bio", k)
        saveRDS(
          readAll(stack(bioclim, r)),
          file = paste0(layer_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds")
        )
      }
    }
  }
}
# unlink(gcm_reg_path, recursive=T)

## Sync NAs
nas <- list()
files <- list.files(layer_path, full.names = TRUE, pattern = region)
stack <- stack(files[grepl(paste0("(?=.*",region,")"), files, perl = TRUE)])
sums <- calc(stack, sum)
sums_nona <- which(is.na(sums[]))

for(i in 1:nlayers(stack)){
  r <- stack[[i]]
  r[sums_nona] <- NA
  writeRaster(r, paste0(layer_path,names(stack)[[i]]), format = "GTiff", overwrite = T)
  print(i)
}
gc()


## Sync NAs (final step for covariates stack and bioq layers) ####
source(file.path(".", "R", "1_functions.R"))

for (region in c("vnm", "aus")){
  layer_path <- file.path(".", "RData")
  
  reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
  covariates <- readRDS(file.path(layer_path, paste0("covariates_", region, ".rds")))
  reg_mask <- covariates[[1]]
  
  files_count  <- list.files(file.path(layer_path), pattern = region, full.names = TRUE)
  bio <- files_count[grepl("bioq", files_count)]
  covs <- files_count[grepl("covariates", files_count)]
  files <- c(bio, covs)
  for(j in 1:length(files)){
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      reg_mask <- mask(reg_mask, r[[i]])
      print(i)
    }
    print(j)
  }
  
  for(j in 1:length(files)){
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      r[[i]] <- mask(r[[i]], reg_mask)
      print(i)
    }
    saveRDS(readAll(r), file = files[[j]])
    print(j)
  }
  print(paste0("processing ", region, "... "))
  length(which(is.na(reg_mask[])))
  length(which(is.na(r[[1]][])))
  which(is.na(reg_mask[]) != is.na(r[[1]][]))
}


## 3. Bioregions layer (for Australia only) ####
## Source: http://www.environment.gov.au/fed/catalog/search/resource/downloadData.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D
bioreg <- readOGR(file.path(regSSP_birds_data, "IBRA7_regions", "ibra7_regions.shp"))
reg_mask <- readRDS(file = file.path(".", "RData", "mask_aus.rds"))
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
saveRDS(bioreg_rast, file.path(".", "RData", "bioregions_aus.rds"))


## 4. GTAP data #### ---- to be updated ----
regions <- c("aus", "vnm")
rcps <- c("1oC", "4oC")
scens <- c("26", "85")
timesteps <- c(1, 13, 23, 33, 43, 53) #(#years from 2018)

## Land endowment (land harvested) for 2016 based on FAO data; available only for 7 commodities(?)
harvested <- read.csv("./raw_data/fao_data/FAOSTAT_data_9-18-2018.csv", header = T)
j <- 2
i <- 2
for (j in 1:length(regions)){
  region <- regions[j]
  landuse <- readRDS(file.path(".", "RData", paste0("covariates_", region, ".rds")))
  landuse <- landuse[[which(names(landuse) == "landuse")]]
  demand_2005 <- table(landuse[])
  
  for (i in 1:length(rcps)){
    rcp <- rcps[i]
    
    if (region == "aus"){
      harv <- harvested[which(harvested$Area == "Australia"),]
    }else if (region == "vnm"){
      harv <- harvested[which(harvested$Area == "Viet Nam"),]
    }
    
    ## sum the area (in ha) produced in available GTAP classes
    total_bysector <- aggregate(harv$Value, by = list(harv$GTAP.sector), FUN = function(x) sum(na.omit(x)))
    
    ## relative area harvested of agricultural commodities as share of total area harvested in one of the GTAP sectors.
    total_bysector[,2] <- total_bysector[,2]/sum(total_bysector[,2])
    saveRDS(total_bysector, file.path(".", "RData", paste0("harvested2016_", region, ".rds")))
    
    
    ## Land endowments for each year from 2019 - 2070 based on GTAP data; availabel for all 58 commodities
    if (region == "aus") {
      cn <- "aus"
    } else if (region == "vnm") {
      cn <- "vn"
    }
    
    c <- read.table(file.path(".", "raw_data", "gtap_data", "rcps",
                              paste0("qfe_", cn, "_", scens[i], "_all.csv")), 
                    header = F, sep = ",")
    c_df <- as.table(t(c[which(c$V1 == "land"),][-1]))
    coms <- t(t(c[2,][1,]))
    colnames(c_df) <- c$V1[seq(1,98, 8)]
    yrs <- c(2019:2070)
    m <- matrix(NA, nrow(c_df), length(yrs))
    for(k in 1:nrow(c_df)){
      ## linear interpolation to estimate land endowments for each year from 2019 - 2070
      m[k,] <- approx(y = c_df[k,-13], x = as.numeric(colnames(c_df))[-13], xout = yrs)$y
    }
    
    m <- as.data.frame(cbind("2018" = 0, m))
    rownames(m) <- coms[-1]
    colnames(m)[-1] <- as.character(yrs)
    ## only keep commodities as per FAO data [7 commodities]
    m_agr <- m[rownames(m)%in%total_bysector$Group.1,]
    m_agr <- m_agr[match(rownames(m_agr), total_bysector[,1]),]
    ## add forestry and cattle
    m_agr <- rbind(m_agr, m[which(rownames(m)%in%c("frs", "ctl")),])
    saveRDS(m_agr, file = file.path(".", "RData", paste0("gtap_landendowments_", region,"_", scens[i], ".rds")))
    ## print land endowments for 2070
    print(m_agr[,ncol(m_agr)])
  }
}

## Difference between total_bysector(fao) & m_agr(gtap)
temp <- t(t(m_agr[,ncol(m_agr)]))
rownames(temp) <- rownames(m_agr)
colnames(temp)[1] <- 'endowment_2070_gtap' #m_agr for 2070
temp <- cbind(temp, c(total_bysector$x, NA, NA))
colnames(temp)[2] <- 'harvested_2016_fao' # total_bysector



## 7. Biodiversity data ####
## Australia:  https://doi.org/10.15468/dl.khlzmu
## Tile 29: https://doi.org/10.15468/dl.yapqxq
## Viet Nam: https://doi.org/10.15468/dl.nt0ftl

## ***** Specify region for analysis (up till Bioclim - future) ***** ####
##  Rerun for c("vnm", "aus", "til")
region <- 'vnm'
if(region == "til"){region <- 'xx'}

source('./R/filter.gbifcsv.R')

## Select country and load GBIF data
gbif_dat <- as.data.frame(fread(paste0(file.path(regSSP_birds_data, "gbif/gbif_aves_"), region, ".csv"), header = T, select = c("decimallongitude", "decimallatitude", "species"), na.strings=c("NA", "", " ")))

gbif_all <- as.data.frame(fread(paste0(file.path(regSSP_birds_data, "gbif/gbif_aves_"), region, ".csv"), header = T, na.strings=c("NA", "", " ")))

gbif_filter <- filter.gbifcsv(gbif_all)

reg_mask <- readRDS(file.path(layer_path, paste0("mask_", region, ".rds")))
## Find exact spatial duplicates and remove
species <- unique(gbif_dat$species)
occ2 <- data.frame()

## Remove dups within same species
for (i in 1:length(species)) {
  spec_inds <- which(occ$species == species[i])
  values <- extract(reg_mask, occ[spec_inds, c(1,2)], cellnumbers = T)
  occ2 <- rbind(occ2, occ[spec_inds[which(!duplicated(values[,1]))],])
  print(i)
}

## Remove entries with less than 2 dec places (because their precision is ca. below 1km2)
lis <- strsplit(sub('0+$', '', as.character(occ2$decimallongitude)), ".", fixed = TRUE)
out <- sapply(lis, function(x) nchar(x[2]))
lis <- strsplit(sub('0+$', '', as.character(occ2$decimallatitude)), ".", fixed = TRUE)
out2 <- sapply(lis, function(x) nchar(x[2]))
unprecise <- which(out < 2 | out2 < 2)

if(length(unprecise) > 0) { occ2 <- occ2[-unprecise,] }

## check that all obs are in locations with data and remove obs where that isn't the case
bgp <- SpatialPoints(occ2[,c(1,2)])
preds_obs <- extract(reg_mask, bgp)
occ2 <- occ2[-which(is.na(preds_obs)),]

## Filter for >= 20 occurrences
occ3 <- occ2[occ2$species%in%names(which(table(occ2$species) >= 20)),]
min(table(occ3$species)) >= 20
species_left <- unique(occ3$species)
spnums <- table(occ3$species)

dat <- occ3
names(dat) <- c("long","lat","species") #column names

## Subset the ones only occuring in vnm
if(region == "til"){
  specs <- unique(dat$species)
  specs <- specs[which(specs%in%unique(occ$species))]
  dat <- dat[dat$species%in%specs,]
}

saveRDS(dat, file = paste0(layer_path, "/occ_", region, ".rds"))


## EXTRAS ###
## Bioclim current - Download and mask
# if(region == "aus"){
#   tiles <- c("39", "310", "311", "49", "410", "411")
# }else if(region%in%c("vnm", "til")){
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
# 
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
