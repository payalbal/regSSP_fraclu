## Processing input data for analysis (using fractional land use)


## Master log file
job_start <- Sys.time()
masterlog <- paste0("./preprocessing_run.txt")
writeLines(c(""), masterlog)
cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)


## Set up work environment ####
setwd("...")

rm(list = ls())
gc()
# devtools::install_github('kapitzas/WorldClimTiles')
# devtools::install_github('skiptoniam/sense')
system(paste0("curl https://raw.githubusercontent.com/skiptoniam/sense/master/R/gdal_raster_functions.R -o ", "./scripts/gdal_raster_functions.R"))
source("./scripts/gdal_raster_functions.R")
# devtools::install_github('smwindecker/gdaltools')
x <- c('data.table','rgdal','rgeos','matrixStats',"sp",'raster',
       'WorldClimTiles','sense' , 'readxl', 'gdaltools', 'sf')
lapply(x, require, character.only = TRUE)
source(file.path(".", "scripts", "0_functions.R"))
gsdms_data <- "/Volumes/uom_data/gsdms_data" #server: "/tempdata/workdir/data"
regSSP_data <- "/Volumes/uom_data/regSSP_data" #server: "/tempdata/workdir/regSSP_data"
rdata_path <- file.path(regSSP_data, "RData") #server: "./RData" 
if(!dir.exists(rdata_path)){dir.create(rdata_path)}
dstrast_path <- file.path(rdata_path, "dst_rasters") #server: "file.path(gsdms_data, "dst_rasters")"
# dir.create(dstrast_path)


## 1. Prepare masks ####
regions <- c("vn", "aus", "til")

for (region in regions){
  
  if (region == "til") {
    reg_mask <- raster(extent(90, 120, 0, 30),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
                       res = 0.008333333)
    reg_mask[] <- 1
    saveRDS(reg_mask, file.path(rdata_path, paste0("mask_", region, ".rds")))
  } else {
    reg_mask <- getData("GADM", country = region, level = 0, path = rdata_path)
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
    saveRDS(reg_mask, file.path(rdata_path, paste0("mask_", region, ".rds")))
    plot(reg_mask)
  }
}
rm(mask_template)


## ***** Specify region for analysis (up till Bioclim - future) ***** ####
##  Rerun for c("vn", "aus", "til")

region <- 'aus'
# for (region in regions){

## 2. Prepare covariate data ####
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))

## Processing at 100km2
reg_mask <- aggregate(reg_mask, fact = 10)


## 2a. Topography ####
## Source: https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
## Select polygon including vn and aus
srtm <- raster(file.path(regSSP_data, "sdat_10008_1_20191110_184812575.asc"))
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

# ## Alternate SRTM source: http://srtm.csi.cgiar.org/srtmdata/ (extent= -180,180,-60,60)
#  ##  on weblink click on link for enntire globe.
#  setwd("../data")
#  system("wget http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_30x30/TIFF/N00E090.zip -O srtm_csi.zip")
#  temp <- unzip("srtm_csi.zip", list = TRUE)
#  unzip("srtm_csi.zip", exdir = "./")
#  file.rename(paste0("./", temp$Name), sub(tools::file_path_sans_ext(basename(temp$Name)),"srtm_csi", paste0("./", temp$Name)))
#  file.remove("srtm_csi.zip")

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
## Rasters showing shortest distance to roads, built-up areas, lakes and rivers

## Roads
## Source: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
# infile <- file.path(regSSP_data, "groads-v1-oceania-east-shp", "gROADS-v1-oceania-east.shp")
# unlink(file.path(dstrast_path, "roads_raster.tif"))
# outfile <- file.path(dstrast_path, "roads_raster.tif")
# gdaltools::rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])
# 
# infile <- file.path(dstrast_path, "roads_raster.tif")
# unlink(file.path(dstrast_path, paste0("dstroads_", region, ".tif")))
# outfile <- file.path(dstrast_path, paste0("dstroads_", region, ".tif"))
# proximity_ras(infile, outfile)
# 
# infile <- file.path(regSSP_data, "groads-v1-oceania-east-shp", "gROADS-v1-oceania-east.shp")
# roads <- sf::st_read(infile)
# road_id <- matrix(c(0:7, "hwy", "pri", "sec", "tert", "loc", "trail", "priv", "unspec"), ncol = 2, nrow = 8)
# road_id <- data.frame(road_id)
# road_classes <- sort(unique(roads$FCLASS))
# road_classes <- road_classes[-c(5, 3)]
# for (i in road_classes){
#   print(paste0("Writing subset ", i))
#   out <- roads[which(roads$FCLASS == i),]
#   outfile <- file.path(dstrast_path, paste0("roads_", road_id[i+1,2], "_raster.shp"))
#   test <- st_crop(out, reg_mask)
#   #ggplot() +  geom_sf(data = boundary, fill = "darkred") + geom_sf(data = test)
#   
#   if(nrow(test) == 0){
#     message("road type not in study area")
#     next}
#   st_write(out, outfile)
#   
#   print(paste0("Rasterizing subset ", i))
#   infile <- outfile
#   outfile <- file.path(dstrast_path, paste0("roads_", road_id[i+1,2], "_raster.tif"))
#   rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])
#   unlink(infile)
#   
#   print(paste0("Calculating subset ", i))
#   infile <- outfile
#   outfile <- file.path(dstrast_path, paste0("dstroads_", road_id[i+1,2], "_" , region, ".tif"))
#   proximity_ras(infile, outfile)
#   unlink(infile)
# }
plot(stack(list.files(dstrast_path, pattern = "dstroads", full.names = TRUE)))


## Built-up areas
## Source: http://ref.data.fao.org/map?entryId=c22837d0-88fd-11da-a88f-000d939bc5d8&tab=metadata
infile <- file.path(gsdms_data, "builtup", "GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0.tif")
infile_res <- crs(raster(infile))
unlink(file.path(dstrast_path, paste0("dstbuiltup_", region, ".tif")))
outfile <- file.path(dstrast_path, paste0("builtup_", region, "_reproj.tif"))
gdalUtils::gdalwarp(infile, outfile, s_srs = infile_res, 
                    t_srs = crs(reg_mask), tr = res(reg_mask))

infile <- outfile
outfile <- file.path(dstrast_path, paste0("builtup_", region, "_crop.tif"))
sense::gdalCrop(inpath = infile, outpath = outfile, extent=extent(reg_mask), 
                resolution=res(reg_mask), return = FALSE)

infile <- outfile
outfile <- file.path(dstrast_path, paste0("dstbuiltup_", region, ".tif"))
sense::gdalDistance(infile, outfile)

  # infile <- outfiile
  # outfile <- file.path(dstrast_path, paste0("dstbuiltup_crop_", region, ".tif"))
  # aus.shp <- rnaturalearth::ne_countries(country = "australia",
  #                                         returnclass = "sp")
  # shapefile(aus.shp, file = file.path(dstrast_path, "aus.shp"))
  # aus.shp <- readOGR(file.path(dstrast_path, "aus.shp"))
  # gdalClipWithShape(inrast = infile, inshp = aus.shp,
  #                   outrast = outfile, cropToShape = TRUE)


## Lakes
## Source: http://www.soest.hawaii.edu/wessel/gshhg/
## See documentation here: /Volumes/uom_data/gsdms_data/lakesrivers/SHAPEFILES.TXT
## Read in files for resolution l, levels 2-4
l2 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS_shp/l", "GSHHS_l_L2.shp"))
  # ## same as
  # l2 <- sf::st_read(l2)
  # l2sp <- as(l2, "Spatial")
l3 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS_shp/l", "GSHHS_l_L3.shp"))
l4 <- readOGR(file.path(gsdms_data, "lakesrivers/GSHHS_shp/l", "GSHHS_l_L4.shp"))

## Combine levels 2-4
temp <- rgeos::gUnion(l2,l3)
temp2 <- rgeos::gUnion(temp,l4)
shapefile(x = temp2, file = file.path(dstrast_path, "lakes_l24.shp"))

## Convert to raster & crop
infile <- file.path(dstrast_path, "lakes_l24.shp")
unlink(file.path(dstrast_path, paste0("lakes_l24_", region, ".tif")))
outfile <- file.path(dstrast_path, paste0("lakes_l24_", region, ".tif"))
rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])
# sense::gdalRasterise(infile, outfile)

## Create proximity raster
infile <- outfile
outfile <- file.path(dstrast_path, paste0("dstlakes_l24_", region, ".tif"))
gdalDistance(infile, outfile)

## Rivers
## Source: http://www.soest.hawaii.edu/wessel/gshhg/
## See documentation here: /Volumes/uom_data/gsdms_data/lakesrivers/SHAPEFILES.TXT
## Read in files for resolution l, levels 2-9
river_files <- list.files(file.path(gsdms_data, "lakesrivers/WDBII_shp/l"), full.names = TRUE)
river_files <- grep(".shp$", river_files, value = TRUE)
river_files <- grep("L02|L03|L04|L05|L06|L07|L08|L09", river_files, value = TRUE)

## Combine levels 2-9
temp <- readOGR(river_files[1])
for(i in 2:length(river_files)) {
    temp <- rgeos::gUnion(temp,readOGR(river_files[i]))
}
shapefile(x = temp2, file = file.path(dstrast_path, "rivers_l29.shp"))

## Convert to raster & crop
infile <- file.path(dstrast_path, "rivers_l29.shp")
unlink(file.path(dstrast_path, paste0("rivers_l29_", region, ".tif")))
outfile <- file.path(dstrast_path, paste0("rivers_l29_", region, ".tif"))
rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])

## Create proximity raster
infile <- outfile
outfile <- file.path(dstrast_path, paste0("dstrivers_l29_", region, ".tif"))
gdalDistance(infile, outfile)


## Protected Areas
pafiles <- list.files(file.path(gsdms_data, "protectedareas/WDPA_Aug2020-shapefile"), pattern = "WDPA_Aug2020-shapefile", full.names = TRUE, recursive = TRUE)
pafiles <- pafiles[grep("polygons.shp$", pafiles)]

## Combine files
temp <- readOGR(pafiles[1])
for(i in 2:length(pafiles)) {
  temp <- rgeos::gUnion(temp,readOGR(pafiles[i]))
}
shapefile(x = temp2, file = file.path(gsdms_data, "processed", "wdpa.shp"))

## Convert to raster & crop
infile <- file.path(gsdms_data, "processed", "wdpa.shp")
unlink(file.path(dstrast_path, paste0("wdpa_", region, ".tif")))
outfile <- file.path(dstrast_path, paste0("wdpa_", region, ".tif"))
rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])

## Create proximity raster
infile <- outfile
outfile <- file.path(dstrast_path, paste0("dstwdpa_", region, ".tif"))
gdalDistance(infile, outfile)



infile <- file.path(raw_path, "Global", "Global Protected Areas", "WDPA_Mar2018-shapefile-polygons.shp")
outfile <- file.path(rdata_path, "PA_cropped.shp")
crop_shp(infile, outfile, ext = extent(reg_mask))
test <- sf::st_read(outfile)
plot(test[which(test$IUCN_CAT%in%c("Ia", "Ib", "II")),], max.plot = 1, add = TRUE)

years <- yrs[ts]
for (i in 1:length(years)){
  PA <- test[which(test$IUCN_CAT%in%c("Ia", "Ib", "II") & test$STATUS_YR <= years[i]),]
  r <- rasterize(PA, reg_mask)
  out <- reg_mask
  out[which(!is.na(r[]))] <- 0
  out <- raster::mask(out, reg_mask)
  plot(out)
  writeRaster(out, file.path(rdata_path, paste0("PA", years[i], "_ama.tif")), format = "GTiff", overwrite = TRUE)
}



## 2e. Population density ####
## Source: http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
popdens <- raster(file.path(gsdms_data, "popdens", "gluds00ag.bil"))
popdens <- crop(popdens, reg_mask, snap = "near")
popdens <- mask(popdens, reg_mask)
names(popdens) <- "popd"


## 2f. Landuse ####
## Source: Copernicus data (for fraclu analysis)
## Link @ GEE: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V_Global
## Direct download link: https://zenodo.org/communities/copernicus-land-cover/search?page=1&size=20
## Data processing for global layer @ GEE by MC Loor: https://code.earthengine.google.com/?scriptPath=users%2Fpayalbal%2Fglobal_layers%3Afraclu_mcloor_sklu. See processing in GEE for reclassification scheme.

## Land use classes used:
## 1	urban
## 2	crop
## 3	forest
## 4	grass
## 5	other 
## 6	NA

# ... gee.py script
# ...download by title

lu <- stack(file.path(regSSP_data, "copernicus_fractional/sklu_classes", paste0("copernicus_resampled_1000_", region, ".tif")))
lu <- projectRaster(lu, reg_mask)
lu <- stack(lu)
lu <- mask(lu, reg_mask)


## 2g. Bioclim - current ####
## Download tiles, mosaic, crop and write
file_in <- list.files(file.path(gsdms_data, 'bio_30s'), full.names = T)
bioclim <- list()

for (f in 1:length(file_in)){
  bioclim[[f]] <- crop(raster(file_in[f]), reg_mask)
  bioclim[[f]] <- mask(bioclim[[f]], reg_mask)
}
bioclim <- stack(bioclim)
names(bioclim) <- paste0("bio", 1:19)  

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
saveRDS(readAll(covariates), file = paste0(rdata_path, "/covariates_", region, ".rds"))
saveRDS(reg_mask, file = paste0(rdata_path, "/mask_", region, ".rds"))

rm(terrain, soil, distances, pa, popdens, bioclim, lu)
rm(reg_mask, covariates, file_in)
gc()
# } ## for region{}

## 2h. Bioclim - future ####
## Global variables
# regions <- c("vn", "aus")
# rcps <- c("45", "60", "85")
# models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")


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
# gcm_reg_path <- file.path(regSSP_data, 'gcm_reg')
gcm_quant_path <- file.path(regSSP_data, 'gcm_quant')
# gcm_reg_path <- file.path(regSSP_data, 'gcm_reg')
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
#         reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
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
#   reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
#   r <- reg_mask
#   inds <- which(!is.na(r[]))
#   
#   j<-1
#   for(j in 1:length(rcps)){
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
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
#       c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
#       for(m in 1:3){
#         bioclim <- readRDS(file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#         r[inds] <- c[,m]
#         names(r) <- paste0("bio", k)
#         saveRDS(readAll(stack(bioclim, r)), file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#       }
#     }
#   }
# }
file.copy(list.files(gcm_quant_path, full.names = TRUE), file.path(rdata_path))
rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
gc()


## Sync NAs - II ### 
## For covariartes + bioq layers
## Find min non-NA set values across mask, covariates and bioq layers and sync NAs 
regions <- c("vn", "aus") ## or c("vn", "aus", "til")
for (region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
  covariates <- readRDS(file.path(rdata_path, paste0("covariates_", region, ".rds")))
  
  files  <- list.files(file.path(rdata_path), pattern = region, full.names = TRUE)
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
  saveRDS(reg_mask, file = paste0(rdata_path, "/mask_", region, ".rds"))
  
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

rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
gc()


# ## Reduce raster size ####
# ##  Remove NA values from covariate stack 
# ##  Mask retained with NA values
# regions <- c("vn", "aus", "til")
# rdata_path2 <- file.path(regSSP_data, "nonatables") # change as per server: boab - "./nonatables"
# if(!dir.exists(rdata_path2)){dir.create(rdata_path2)}
# 
# for(region in regions){
#   print(paste0("processing ", region, "... "))
#   reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
#   reg_layers <- list.files(rdata_path, full.names = TRUE, pattern = region)
#   ## files[grepl(paste0("(?=.*",region,")"), files, perl = TRUE)]
#   reg_layers <- reg_layers[!grepl("mask", reg_layers)]
#   
#   for(j in 1:length(reg_layers)){
#     print(paste0("processing file = ", j, " of ", length(reg_layers) ,":", reg_layers[j]))
#     reg_stack <- readRDS(reg_layers[j])
#     ind_nona <- which(!is.na(reg_mask[]))
#     
#     print(paste0("# layers in ", reg_layers[j], " = ", nlayers(reg_stack)))
#     new_dat <- data.table()
#     for(i in 1:nlayers(reg_stack)){
#       print(paste0("processing layer = ", i, " of ", nlayers(reg_stack)))
#       r <- reg_stack[[i]]
#       r <- getValues(r)
#       r <- r[ind_nona]
#       new_dat[,paste0("new_dat", "$", names(reg_stack[[i]])) := r]    
#     }
#     saveRDS(new_dat, file = paste0(rdata_path2, basename(reg_layers[j]), "_nona.", "rds"))
#     ## txt or csv files using fwrite much larger!
#   }
# }
# 
# rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
# gc()

## Save raster stacks as data.table ####



## 3. Bioregions layer (for Australia only) ####
## Source: http://www.environment.gov.au/fed/catalog/search/resource/downloadData.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D
bioreg <- readOGR(file.path(regSSP_data, "IBRA7_regions", "ibra7_regions.shp"))
reg_mask <- readRDS(file = file.path(rdata_path, "mask_aus.rds"))
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
saveRDS(bioreg_rast, file.path(rdata_path, "bioregions_aus.rds"))


## 4. GTAP data ####
regions <- c("aus", "vn")
ssps <- paste0("ssp", 1:3)

## Area harvested for crops in 2016 - FAO data
## Data source:http://www.fao.org/faostat/en/#data/QC
## 'Item' in data is coded for GTAP sectors/commodities befor importing in R
harvested <- read.csv(file.path(regSSP_data,"fao_data/FAOSTAT_data_9-18-2018.csv"), header = T) 

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
  saveRDS(total_bysector, file.path(rdata_path, paste0("harvested2016_", region, ".rds")))
  
  ## Specify additional gtap sectors (i.e. forestry, cattle)
  gtap_sectors <- c(as.character(total_bysector$gtap_sector), "frs", "ctl")
  
  ## Load GTAP land endowmenet by sector data for 2020 - 2071
  ## Source: Tom Kompas (Uni Melb), Van Ha Pham (ANU)
  for (ssp in ssps){
    temp <- as.data.table(read_xls(file.path(regSSP_data, "gtap_data", "ssps",
                                             paste0(region, "_", ssp, ".xls")), sheet = "land"))
    # temp <- temp[X__1 %in% gtap_sectors]
    temp <- temp[...1 %in% gtap_sectors] # on server X__1 is renamed as ...1
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
    saveRDS(new_temp, file = file.path(rdata_path, paste0("gtap_landendowments_", region,"_", ssp, ".rds")))
    
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
source(file.path(".", "scripts", "0_filter.gbifcsv.R"))
gbif_backbone <- fread(file.path(regSSP_data, "gbif","Taxon.tsv"))
regions = c("aus", "til")
# For vn, we fit models to larger data from til (subset by vnm species) 
vnm_species <- unique(fread(file.path(regSSP_data, "gbif", 
                                      paste0("gbif_aves_vn.csv")), header = T, 
                            na.strings=c("NA", "", " "))$species)

for(region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
  
  gbif_all <- fread(file.path(regSSP_data, 
                              "gbif", paste0("gbif_aves_", region, ".csv")), 
                    header = T, na.strings=c("NA", "", " "))
  
  filter.gbifcsv(gbif.downloaded.data = gbif_all,
                 gbif.nub.taxonomy = gbif_backbone,
                 subset.gbifnubtaxonomy.byclass = "Aves",
                 output_folder = rdata_path,
                 output_name = paste0("occ_", region, ".rds"),
                 domain.mask = reg_mask,
                 select_fields = c("decimallongitude", "decimallatitude", "species"))
  
  print(paste0("file saved for ", region, ": ", file.path(rdata_path, paste0("occ_", region, ".rds"))))
  
  ## Subset species for tile 29 as per species in Viet Nam
  if(region == "til"){
    dat_til <- readRDS(file.path(rdata_path, paste0("occ_til.rds")))
    til_species <- unique(dat_til$species)
    til_species <- til_species[which(til_species%in%vnm_species)]
    dat_til <- dat_til[dat_til$species%in%til_species,]
    saveRDS(dat_til, file = file.path(rdata_path, paste0("occ_til.rds")))
  }
}


## Log job duration
job_end <- Sys.time()
job_duration = job_end - job_start
cat(paste0(">> Job end = ", job_end, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Job duration for region ", region, " = ", job_duration, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Date:\n"), file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
Sys.Date()
sink(NULL)
cat("\n\nSession Info:\n", file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
sessionInfo()
sink(NULL)

rm(list=ls())
gc()


