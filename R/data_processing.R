
## Data download and processing for regional trade analysis (Vnm + Aus)

## Set work environment
x <- c('data.table', 'rgdal', 'rgeos', 'sp', 'raster', 'rgbif', 'sp', 'GADMTools', 'sf', 'tools')
lapply(x, require, character.only = TRUE)
# devtools::install_github('skiptoniam/sense')
library(sense) 


## Create folders
# if(!dir.exists("data")){dir.create("data")} # on server
# data_path <- "./data" # on server

data_path <- "/Volumes/discovery_data/regtrade_data" # on local machine
if(!dir.exists(file.path(data_path, "processed"))){dir.create(file.path(data_path, "processed"))}


## GBIF backbone taxonomy
##  GBIF Secretariat (2017). GBIF Backbone Taxonomy. Checklist dataset https://doi.org/10.15468/39omei accessed via GBIF.org on 2019-08-26.
system(paste0("curl http://rs.gbif.org/datasets/backbone/backbone-current.zip -o ", data_path, "/gbif_taxonomy.zip"))
unzip(file.path(data_path, "gbif_taxonomy.zip"), list = TRUE)
unzip(file.path(data_path, "gbif_taxonomy.zip"), files = 'Taxon.tsv', exdir = data_path)
backbone <-fread("data/Taxon.tsv")


## GBIF occurrence data
get_gbif <- function (folder_path, url, file_name) {
  system(paste0("curl ", url, " -o ", file.path(folder_path, file_name)))
  temp <- unzip(file.path(folder_path, file_name), list = TRUE)
  unzip(file.path(folder_path, file_name), files = 'occurrence.txt', 
        exdir = folder_path)
  file.rename(file.path(data_path, "occurrence.txt"), 
              paste0(data_path, "/", file_path_sans_ext(file_name), ".txt"))
  file.remove(file.path(folder_path, file_name))
}

##  Vietnam birds
##  In GBIF select Aves in 'Scientific Name' and Viet Nam in 'Country of area'
##  GBIF.org (30 August 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.sjdkid
url = "http://api.gbif.org/v1/occurrence/download/request/0007376-190813142620410.zip"
file_name = "vnm_gbif.zip"
get_gbif (data_path, url, file_name)

system(paste0("curl http://api.gbif.org/v1/occurrence/download/request/0007376-190813142620410.zip -o ", data_path, "/vnm_gbif.zip"))
temp <- unzip(file.path(data_path, "vnm_gbif.zip"), list = TRUE)
unzip(file.path(data_path, "vnm_gbif.zip"), files = "occurrence.txt", exdir = data_path)
file.rename(file.path(data_path, "occurrence.txt"), "vnm_gbif.txt"))
file.remove(file.path(data_path, "vnm_gbif.zip"))


##  Tile29 birds
##  In GBIF select Aves in 'Scientific Name' and specfy 'Location' as 
##  85, 125, -5, 35 (xmin, xmax, ymin, ymax) i.e. extent of til29_mask + 5
##  GBIF.org (30 August 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.sjdkid
system(paste0("curl http://api.gbif.org/v1/occurrence/download/request/0008220-190813142620410.zip -o ", data_path, "/til29_gbif.zip"))
temp <- unzip(file.path(data_path, "til29_gbif.zip"), list = TRUE)
system(paste0("unzip ", file.path(data_path, "til29_gbif.zip"), " -d ", data_path)) ## 4GB limit on unzip within R
file.rename(file.path(data_path, temp$Name), sub(file_path_sans_ext(basename(temp$Name)),"til29_gbif", file.path(data_path, temp$Name)))
file.remove(file.path(data_path, "til29_gbif.zip"))

  ## Vietnam mammals: http://api.gbif.org/v1/occurrence/download/request/0007378-190813142620410.zip                              
  ## GBIF.org (30 August 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.ncsbla
  
  ## Australia birds: http://api.gbif.org/v1/occurrence/download/request/0007392-190813142620410.zip
  ## GBIF.org (30 August 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.4agmjw
  
  ## Australia mammals: http://api.gbif.org/v1/occurrence/download/request/0007396-190813142620410.zip
  ## GBIF.org (30 August 2019) GBIF Occurrence Download https://doi.org/10.15468/dl.vfdbud


## Covariate data
##  SRTM
##  Source: http://srtm.csi.cgiar.org/srtmdata/
##  on weblink select 30 x 30 tile size and Geo TIFF format > Search
system("wget http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_30x30/TIFF/N00E090.zip -O data/srtm.zip")
temp <- unzip("data/srtm.zip", list = TRUE)
unzip("data/srtm.zip", exdir = "./data")
file.rename(paste0("./data/", temp$Name), sub(file_path_sans_ext(basename(temp$Name)),"srtm", paste0("./data/", temp$Name)))
file.remove("data/srtm.zip")

##  Soil
##  Source (needs login): https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
##  on weblink check the files you want to download and place an order, then copy link from email recieved
system("wget -r -np -nH --reject 'index.html*' -e robots=off https://daac.ornl.gov/orders/ca2c9c7d67f425992259251a578b8669/ -O data/soil.zip")
temp <- unzip("data/soil.zip", list = TRUE)
unzip("data/soil.zip", exdir = "./data")
file.rename(paste0("./data/", temp$Name), sub(file_path_sans_ext(basename(temp$Name)),"soil", paste0("./data/", temp$Name)))
file.remove("data/soil.zip")

## Roads
## Source: 
https://sedac.ciesin.columbia.edu/downloads/data/groads/groads-global-roads-open-access-v1/groads-v1-global-gdb.zip


area_sdm <- c("til29", "vnm")

for (area_sdm in area_sdm){
  
  ## Load mask and gbif data
  if (area_sdm == "til29") {
    reg_mask <- getData('worldclim', var = 'alt', res=0.5, lon=95, lat=15, path = "./data")
    reg_mask[which(!is.na(reg_mask[]))] <- 1
    gbif.raw <- fread(file.path("data/", "til29_gbif.csv"), header = T, na.strings=c("NA", "", " "))
  } else {
    reg_mask <- getData('alt', country = 'VNM', mask = TRUE, path = "./data")
    reg_mask[which(!is.na(reg_mask[]))] <- 1
    gbif.raw <- fread(file.path("data/", "vnm_gbif.csv"), header = T, na.strings=c("NA", "", " "))
  }
  
    ## Filter GBIF data and save as reduced data.table
  source("./R/00_filter.gbifraw.R")
  filter_gbifraw(gbif.downloaded.data = gbif.raw, gbif.nub.taxonomy = backbone, 
                    output_folder = "data/processed", output_name = paste0("filtered_gbif_", area_sdm), 
                    domain.mask = reg_mask,
                    start.year = 1950, end.year = 2018,
                    spatial.uncertainty.m = 1000,
                    verbose = TRUE)
  
  rm(gbif.raw)
  gc()
  

  ## Load predictor data
  ## SRTM Global 30 arc secs NASA;	https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
  srtm <- raster(file.path("data/raw", "srtm_aus_vnm_tile29.tif"))
  srtm <- crop(srtm, reg_mask)
  names(srtm) <- "srtm"
  extent(srtm) <- extent(reg_mask)
  elevation <- mask(srtm, reg_mask)
  slope <- terrain(srtm, opt = "slope")
  roughness <- terrain(srtm, opt = "roughness")
  aspect <- terrain(srtm, opt = "aspect")
  terrain <- stack(elevation, slope, roughness, aspect)
  
  for (i in 1:3){
    r <- terrain[[i]]
    r <- projectRaster(r, reg_mask)
    r <- mask(r, reg_mask)
    print(cellStats(r, mean))
    writeRaster(r, paste0(data_path, "/", paste(c(names(terrain)[i], area_sdm, "sta","000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  }
  
  rm(r, srtm, elevation, slope, roughness, aspect)
  gc()
  
  ## Soil varaibales
  ## source: Global Soil Data Task Group;	https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
  names <- c( "soilbulkdens", "soilavailablewataer", "soilcarbondens",  "soilnitro")
  soil <- stack(list.files(file.path("data/raw", "soil", "data"), pattern = "*.dat", full.names = T))
  
  for (i in 1:4){
    r <- raster(soil[[i]])
    crs(r) <- crs(reg_mask)
    r <- projectRaster(r, reg_mask)
    r <- mask(r, reg_mask)
    writeRaster(r, paste0(data_path, "/", paste(c(names[i], area_sdm, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  
  ## 'Distance to' variables
  ## sources: 
  ## Protected areas: https://www.protectedplanet.net/c/about
  ## Road network: Socioecnomic data and application centre; http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
  ## Built-up: 
  ## Drainage, Sea, Lakes: NOAA - Global Self-consistent Hierarchical High-resolution Geography;	http://www.soest.hawaii.edu/wessel/gshhg/
  
  pa <- shapefile("data/raw/protectedarea/WDPA_Nov2018-shapefile-polygons.shp")
  pa <- crop(pa, reg_mask) # to clip vector to ratser extent
  r <- raster(extent(reg_mask))
  crs(r) <- crs(reg_mask)
  r <- projectRaster(r, reg_mask, method = "ngb")
  dd = gDistance(pa, as(r,"SpatialPoints"), byid=TRUE)
  r[] = apply(dd,1,min)
  writeRaster(r, paste0(data_path, "_new", "/", paste(c("distpas", area_sdm, "sta", "000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  # pas <- raster("data/processed/layers_sdm_new/distpas_til_sta_000_00.tif") ## does not work, use the code below instead
  dir <- paste0(data_path, "_new", "/")
  files<-list.files(dir, pattern = ".tif")
  pas <- lapply(paste0(dir, files), raster)
  plot (pas[[1]])
  
  
  
  
  distances <- stack(list.files(file.path("data/raw", "qgis_preprocessed"),  full.names = T, pattern = area_sdm))
  distances[is.na(reg_mask[])] <- NA
  ## OR? distances <- mask(distances, reg_mask)??
  names <- c("distbuiltup", "distcoastline", "distlakes", "distrivers", "distroads", "distpas")
  ## ?? pa[] <- (pa[] * -1) + 1
  for(i in 1:nlayers(stack)){
    writeRaster(stack[[i]], paste0(data_path, "/",paste(c(names[i], area_sdm, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  ## --------------------------------------------------------------------------------------------------------------------------------------------
  ## steps for pre-processing:
  ## combine shapefiles (different levels) in qgis (Levels 2-9 for rivers)
  ## rasterize in qgis (automatically subsetted from global data set), because r is too slow for that
  ## proximity in qqis (r too slow)
  ## then load here to mask NAs, align with other layers and write
  
  ## this is so qgis can read the attributes for some of the croppsed shapefiles
  ## have done this already, only needs done once, so code just for documentation purposes
  
  # name <- "WDBII_rivers_global_L2-L9"
  # area_sdm <- "vnm"
  # shp <- readOGR("/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name)
  # shp@data[,2] <- as.numeric(shp@data[,2])
  # shp@data[,1] <- as.numeric(shp@data[,1])
  # writeOGR(shp, dsn = "/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name, driver = "ESRI Shapefile", overwrite = T)
  ## load qgis processed files, and set NA
  # lakes <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/lakes_", area_sdm, ".tif"))
  # coast <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/coast_", area_sdm, ".tif"))
  # rivers <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/rivers_", area_sdm, ".tif"))
  # PA <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", area_sdm, ".tif"))
  # roads <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", area_sdm, ".tif"))
  # builtup <- raster(paste0("/Users/simon/Dropbox/PhD/Chapter 1/processed data/raw_rasters_rivers_lakes_cl/PA_", area_sdm, ".tif"))
  # names <- c("distrivers", "distlakes", "distcoastline", "PA", "distbuiltup", "distroads")
  ## --------------------------------------------------------------------------------------------------------------------------------------------
  
  
  ## Population density
  ## source: Global Rural-Urban Mapping Project, Version 1; http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
  r <- raster(file.path("data/raw", "Pop_density", "gluds00ag.bil"))
  pop_dens <- crop(r, reg_mask, snap = "near")
  extent(pop_dens) <- extent(reg_mask)
  pop_dens <- mask(pop_dens, reg_mask)
  writeRaster(pop_dens, paste0(data_path, "/", paste(c("popdens", area_sdm, "sta", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  rm(r).
  gc()
  
  ## Landuse 
  ## source: Global Land Cover Characterization (GLCC); https://www.usgs.gov/media/images/global-land-cover-characteristics-data-base-version-20 
  world_lu <- raster(file.path("data/raw", "LCType.tif"))
  lu <- crop(world_lu, reg_mask)
  lu <- projectRaster(lu, reg_mask, method = "ngb")
  lu <- mask(lu, reg_mask)
  rm(world_lu)
  gc()
  
  ## Reclassify land use classes
  ## Original land use classes in GLCC data: 
  ##    # 0	Water
  ##    # 1	Evergreen Needle leaf Forest
  ##    # 2	Evergreen Broadleaf Forest
  ##    # 3	Deciduous Needle leaf Forest
  ##    # 4	Deciduous Broadleaf Forest
  ##    # 5	Mixed Forests
  ##    # 6	Closed Shrublands
  ##    # 7	Open Shrublands
  ##    # 8	Woody Savannas
  ##    # 9	Savannas
  ##    # 10	Grasslands
  ##    # 11	Permanent Wetland
  ##    # 12	Croplands
  ##    # 13	Urban and Built-Up
  ##    # 14	Cropland/Natural Vegetation Mosaic
  ##    # 15	Snow and Ice
  ##    # 16	Barren or Sparsely Vegetated
  
  lu[lu[]%in%c(12,14)] <- 100       # 0 crop
  lu[lu[]%in%c(6,7,8,9,10)] <- 200  # 1 grass/shrub
  lu[lu[]%in%c(1,2,3,4,5)] <- 300   # 2 forest
  lu[lu[]%in%c(13)] <- 400          # 3 urban/artificial
  lu[lu[]%in%c(11,15,16,0)] <- 500  # 4 barren
  lu <- lu/100 - 1
  table(lu[])
  names(lu) <- "landuse"
  writeRaster(lu, paste0(data_path, "/", paste(c("landuse", area_sdm, "luc", "000", "00"), collapse = "_"), ".tif"), format = "GTiff", overwrite = T)
  
  
  ## DYNAMIC PREDICTOR DATA
  ## ------------------------------
  ## Bioclim data - present 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  
  ## SK's code to downlaod and merge worldclim data by tile
  # source("./R/00_functions.R")
  # tiles <- "29"
  # temp_folder <- file.path(".", "temp")
  # dir.create(temp_folder)
  # bioclim <- get_wctiles(tile = tiles, var = "bio", path = temp_folder)
  # bioclim <- merge_wctiles(bioclim)
  # names(bioclim) <- paste0("bio", c(1:19))
  # bioclim <- crop(bioclim, reg_mask)
  # bioclim <- mask(bioclim, reg_mask)
  # unlink(temp_folder, recursive = T)
  # covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)
  # saveRDS(covariates, file = paste0(data_path, "covariates_", area_sdm, ".rds"))
  
  
  bionames <- gsub('.{7}$', '', list.files(file.path("data/raw", "bioclim", "present_tile29"), pattern = ".bil"))
  tile29_bio <- stack(list.files(file.path("data/raw/", "bioclim", "present_tile29"), pattern = ".bil", full.names = T))
  crs(tile29_bio) <- crs(reg_mask)
  files_sub <- crop(tile29_bio, reg_mask)
  files_sub <- mask(files_sub, reg_mask)
  for(i in 1:nlayers(files_sub)){
    writeRaster(files_sub[[i]], paste0(data_path, "/", paste(c(bionames[i], area_sdm, "dyn", "000", "00"), collapse = "_"), ".tif"), overwrite = T)
  }
  
  rcp <- c(26, 85)
  ## Bioclim data - 2070, RCP 2.6 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  files <- stack(list.files(file.path("data/raw", "bioclim", "he26bi70"), pattern = "tif", full.names = T, recursive = T))
  files <- crop(files, reg_mask)
  files <- mask(files, reg_mask)
  for(i in 1:nlayers(files)){
    writeRaster(files[[i]], paste0(data_path, "/", paste(c(bionames[i], "26",area_sdm, "dyn", "000", "56"), collapse = "_"), ".tif"), overwrite = T, format = "GTiff")
  }
  
  ## Bioclim data - 2070, RCP 8.5 
  ## source: downloaded for tile29; http://www.worldclim.org/cmip5_30s
  files <- stack(list.files(file.path("data/raw", "bioclim", "he85bi70"), pattern = "tif", full.names = T, recursive = T))
  files <- crop(files, reg_mask)
  files <- mask(files, reg_mask)
  for(i in 1:nlayers(files)){
    writeRaster(files[[i]], paste0(data_path, "/", paste(c(bionames[i], "85", area_sdm, "dyn", "000", "56"), collapse = "_"), ".tif"), overwrite = T, format = "GTiff")
  }
  
  
  
  ## SYNC NAs FOR REDICTOR LAYERS
  ## ------------------------------
  nas <- list()
  files <- list.files(data_path, full.names = T, pattern = area_sdm)
  pred.stack <- stack(files[grepl(paste0("(?=.*",area_sdm,")"), files, perl = TRUE)])
  sums <- calc(pred.stack, sum)
  sums_nona <- which(is.na(sums[]))
  
  for(i in 1:nlayers(pred.stack)){
    r <- pred.stack[[i]]
    r[sums_nona] <- NA
    writeRaster(r, paste0(data_path, "/", names(pred.stack)[[i]]), format = "GTiff", overwrite = T)
    print(paste0(i, "/", nlayers(pred.stack)))
  }
  rm(r, pred.stack, sums, sums_nona)
  gc()
  
  
  
  ## CREATE PREDICTOR LAYERS FOR ALL TIME STEPS VIA LINEAR INTERPOLATION
  #-------------------------------------------------------------------------
  
  ## Load present climate data
  pres_clim <- stack(list.files(data_path, pattern = paste(c( area_sdm,"dyn", "000", "00"), collapse = "_"), full.names = T))
  pres_clim_df <- as.data.frame(pres_clim)
  time_steps <- c(16, 26, 36, 46) # from 2015 onwards: 2030, 2040, 2050, 2060
  ras_template <- pres_clim[[1]]
  
  ## Loop through scenarios and produce linearly interpolated covariates for 10 yr time steps
  for (i in 1:length(rcp)){
    fut_clim <- stack(list.files(data_path, pattern = paste(c(rcp[i], area_sdm, "dyn", "000", "56"), collapse = "_"), full.names = T))
    fut_clim_df <- as.data.frame(fut_clim)
    incr_df <- (fut_clim_df - pres_clim_df)/ (2070 - 2015)
    for (j in 1:length(time_steps)){
      ts_clim_df <- pres_clim_df + time_steps[j] * incr_df
      for (k in 1:nlayers(pres_clim)){
        ras_template[] <-  ts_clim_df[,k]
        writeRaster(ras_template, paste0(data_path,"/", paste(c(bionames[k], "_", rcp[i], area_sdm, "dyn", "000", time_steps[j]), collapse = "_")), format = "GTiff", overwrite = T)
      }
    } # end time_steps
  } # end rcp
} # end area_sdm

rm(backbone)
gc()
