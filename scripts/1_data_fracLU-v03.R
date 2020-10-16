## Processing input data for land use analysis (using fractional land use)
##
## NOTES:
## See gsdms.Rproj for data download scripts
##
## TO FIX:
## Specify correct crs for vietnam...
## Specify pixel size for 10km resolution for vn...


## Master log file ####
job_start <- Sys.time()
masterlog <- paste0("./preprocessing_run.txt")
writeLines(c(""), masterlog)
cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)



## Set up work environment ####
rm(list = ls())
gc()

setwd("./regSSP_fraclu/")

# devtools::install_github('smwindecker/gdaltools')
# devtools::install_github('skiptoniam/sense')
system(paste0("curl https://raw.githubusercontent.com/skiptoniam/sense/master/R/gdal_raster_functions.R -o ",
              "./scripts/gdal_raster_functions.R"))
source("./scripts/gdal_raster_functions.R")

# ## jgarber's repo:
# library(RCurl)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R"
# url.exists(url)
# source(url)
# devtools::source_url(url)
# url <- "https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons.git"
# devtools::install_git(url = url)
source("./gdal_calc.R") # by jgarber

source(file.path(".", "scripts", "0_functions.R"))

x <- c('data.table','rgdal','rgeos','matrixStats',"sp",'raster',
       'readxl', 'gdaltools', 'sf')
lapply(x, require, character.only = TRUE)
rm(x)

# ## Server paths
# setwd("/tempdata/workdir/regSSP_fraclu/")
# gsdms_data <- "/tempdata/workdir/data"
# regSSP_data <- "/tempdata/workdir/regSSP_data"
# rdata_path <- "./RData" 
# if(!dir.exists(rdata_path)){dir.create(rdata_path)}
# dstrast_path <- file.path(gsdms_data, "dst_raster")
# # dir.create(dstrast_path)


## Local paths
gsdms_data <- "/Volumes/uom_data/gsdms_data"
regSSP_data <- "/Volumes/uom_data/regSSP_data"
rdata_path <- file.path(regSSP_data, "RData")
if(!dir.exists(rdata_path)){dir.create(rdata_path)}
dstrast_path <- file.path(rdata_path, "dst_rasters")
if(!dir.exists(dstrast_path)){dir.create(dstrast_path)}



## >> Prepare masks: In WGS84 ####
regions <- c("vn", "aus", "til")

for (region in regions){
  
  if (region == "til") {
    reg_mask <- raster(extent(90, 120, 0, 30),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                       res = 0.008333333)
    reg_mask[] <- 1
    writeRaster(reg_mask, file.path(rdata_path, paste0("mask_", region, ".tif")))
  } 
  
  else {
    reg_mask <- getData("GADM", country = region, level = 0, path = rdata_path)
    reg_mask <- gSimplify(reg_mask, tol = 0.008333333)
    
    if (region == "vn") {
      mask_template <- raster(
        nrow = 1779,
        ncol = 879,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
        ext = extent(c(102.1417, 109.4667, 8.566667, 23.39167))
      )
      mask_template[] <- 1
    }
    
    else if (region == "aus") {
      mask_template <- raster(
        nrow = 4091,
        ncol = 4990,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
        ext = extent(c(112.4667, 154.05,-44.04167,-9.95))
      )
      mask_template[] <- 1
    }
    reg_mask <- mask(mask_template, reg_mask)
    writeRaster(reg_mask, file.path(rdata_path, paste0("mask_", region, ".tif")))
    plot(reg_mask)
  }
}
rm(mask_template)



## Specify region for analysis ####
##  Rerun for c("vn", "aus", "til")
region <- 'aus'


## Load mask file ####
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))

# ## Specify desired projection and resolution for analyses ####
# infile <- file.path(rdata_path, paste0("maskWGS_", region, ".tif"))
# old_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# 
# if (region == "aus"){
#   new_crs <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
#   # reso = c(1000, 1000) # at 1km
#   reso = c(10000, 10000) # at 10km
# } else {
#   new_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#   res0 = c(0.008333333,0.008333333) # at 1km
# }
# 
# outfile <- sub("WGS","", infile)
# gdalwarp(infile, outfile, s_srs = old_crs, t_srs = old_crs, tr = reso) 
# raster(infile)
# raster(outfile)



## >> Prepare proximity rasters ####
## Rasters showing shortest distance to features

## Roads ####
## Source: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
infile <- file.path(regSSP_data, "groads-v1-oceania-east-shp", "gROADS-v1-oceania-east.shp")
unlink(file.path(dstrast_path, "roads_raster.tif"))
outfile <- file.path(dstrast_path, "roads_raster.tif")
gdaltools::rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])

infile <- file.path(dstrast_path, "roads_raster.tif")
unlink(file.path(dstrast_path, paste0("dstroads_", region, ".tif")))
outfile <- file.path(dstrast_path, paste0("dstroads_", region, ".tif"))
gdaltools::proximity_ras(infile, outfile)

dstroad <- outfile
# list.files(dstrast_path, pattern = paste0("dstroads_", region), full.names = TRUE)

# ## Roads raster by road type
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
# dstroad <- list.files(dstrast_path, pattern = "dstroads", full.names = TRUE)


## Built-up areas ####
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

dstbuiltup <- outfile

# infile <- outfiile
# outfile <- file.path(dstrast_path, paste0("dstbuiltup_crop_", region, ".tif"))
# aus.shp <- rnaturalearth::ne_countries(country = "australia",
#                                         returnclass = "sp")
# shapefile(aus.shp, file = file.path(dstrast_path, "aus.shp"))
# aus.shp <- readOGR(file.path(dstrast_path, "aus.shp"))
# gdalClipWithShape(inrast = infile, inshp = aus.shp,
#                   outrast = outfile, cropToShape = TRUE)


## Lakes ####
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

dstlakes <- outfile


## Rivers ####
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

dstrivers <- outfile


## Protected Areas ####
pafiles <- list.files(paste0(gsdms_data, "/protectedareas"), 
                      pattern = "WDPA_Aug2020-shapefile", full.names = TRUE, recursive = TRUE)
pafiles <- pafiles[grep("polygons.shp$", pafiles)]

## Crop
for (i in 1:length(pafiles)){
  infile <- pafiles[i]
  outfile <- paste0(dstrast_path, "/wdpa_cropped_", i, ".shp")
  sense::ogrBBoxClip(infile, outfile, ext = raster::extent(reg_mask), returnShape = FALSE)
}

test <- rgdal::readOGR(outfile[1])
plot(test)

## Subset WDPA layer by IUCN categories and save by year
pa.iucn.cat <- c("Ia", "Ib", "II")
pa.status.yrs <- c("2005", "2010", "2018")

infile <- list.files(dstrast_path, pattern = "wdpa_cropped", full.names = TRUE, recursive = TRUE)
infile <- infile[grep(".shp$", infile)]
outfile <- gsub("_cropped", "", infile)

for (i in 1:length(infile)){
  test <- sf::st_read(infile[i])
  # plot(test[which(test$IUCN_CAT %in% pa.iucn.cat),], max.plot = 1)
  
  for (j in 1:length(pa.status.yrs)){
    temp.outfile <- gsub(".shp$", paste0("_subset_yr", pa.status.yrs[j], "_", region, ".shp"), outfile[i])
    pa <- test[which(test$IUCN_CAT %in% pa.iucn.cat & test$STATUS_YR <= pa.status.yrs[j]),]
    sf::st_write(pa, dsn = temp.outfile, overwrite = TRUE)
  }
}

## Convert to raster
infile <- list.files(dstrast_path, pattern = "_subset_", full.names = TRUE, recursive = TRUE)
infile <- infile[grep(".shp$", infile)]
outfile <- gsub(".shp$", ".tif", infile)
# gdaltools::rasterize_shp(infile, outfile, res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)])
mapply(gdaltools::rasterize_shp, input_file = infile, output_file = outfile, 
       MoreArgs = list(res = res(reg_mask)[1], ext = extent(reg_mask)[c(1,2,3,4)]))

## Combine rasters
infile <-  outfile

for (i in 1:length(pa.status.yrs)) {
  message(cat("Processing year..", pa.status.yrs[i]))
  temp.infile <- grep(pa.status.yrs[i], infile, value = TRUE)
  outfile <- gsub("_1", "", temp.infile[1])
  
  ## Make a template raster file to build onto. Think of this a big blank canvas to add tiles to.
  e <- extent(reg_mask)
  template <- raster(e)
  projection(template) <- crs(reg_mask)
  unlink(outfile)
  writeRaster(template, file = outfile, format="GTiff")
  
  ## Merge all raster tiles into one big raster.
  gdalUtils::mosaic_rasters(gdalfile = temp.infile, 
                            dst_dataset = outfile, of="GTiff")
  gdalinfo(outfile)
}

## Create proximity raster
infile <- list.files(dstrast_path, pattern = "wdpa_subset_", 
                     full.names = TRUE, recursive = TRUE)
outfile <- file.path(dstrast_path, 
                     paste0(gsub("wdpa_subset", "dstwdpa", 
                                 tools::file_path_sans_ext(basename(infile))), 
                            ".tif"))
mapply(sense::gdalDistance, infile, outfile)

dstwdpa <- outfile

## List proximity raster files ####
prox.files <- c(dstroad, dstbuiltup, dstlakes, dstrivers, dstwdpa)

rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data", "reg")))
gc()



## >> Prepare covariate data ####

## Topography ####
## Source: https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
## Select polygon including vn and aus
srtm <- file.path(regSSP_data, "sdat_10008_1_20191110_184812575.asc")

  # ## Alternate SRTM source: http://srtm.csi.cgiar.org/srtmdata/ (extent= -180,180,-60,60)
  #  ##  on weblink click on link for enntire globe.
  #  setwd("../data")
  #  system("wget http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_30x30/TIFF/N00E090.zip -O srtm_csi.zip")
  #  temp <- unzip("srtm_csi.zip", list = TRUE)
  #  unzip("srtm_csi.zip", exdir = "./")
  #  file.rename(paste0("./", temp$Name), sub(tools::file_path_sans_ext(basename(temp$Name)),"srtm_csi", paste0("./", temp$Name)))
  #  file.remove("srtm_csi.zip")


## Soil ####
## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
soil <- list.files(file.path(gsdms_data, "orders"), 
                   pattern = "*.dat", full.names = T, recursive = T)


## Population density ####
## Source: http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
popdens <- file.path(gsdms_data, "popdens", "gluds00ag.bil")


## Landuse ####
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
lu1 <- file.path(file.path(gsdms_data, "copernicus", "global_frac/sklu_classes", "landuse_class1.tif"))
lu2 <- file.path(file.path(gsdms_data, "copernicus", "global_frac/sklu_classes", "landuse_class2.tif"))
lu3 <- file.path(file.path(gsdms_data, "copernicus", "global_frac/sklu_classes", "landuse_class3.tif"))
lu4 <- file.path(file.path(gsdms_data, "copernicus", "global_frac/sklu_classes", "landuse_class4.tif"))
lu5 <- file.path(file.path(gsdms_data, "copernicus", "global_frac/sklu_classes", "landuse_class5.tif"))
landuse <- c(lu1, lu2, lu3, lu4, lu5)


## Bioclim - current ####
bio_current <- list.files(paste0(gsdms_data, "/bio_30s"), 
                          pattern = "bio_current*", full.names = TRUE)
# bio_names <- gsub("_current_", "", bio_current)


## Bioclim - future ####
# regions <- c("vn", "aus")
# rcps <- c("45", "60", "85")
models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")
gcms <- list.files(file.path(gsdms_data, "gcm_30s"), full.names = TRUE)

  # ## Files by RCP
  # bio_rcp45 <- list.files(file.path(gsdms_data, "gcm_30s"), pattern = "*bc45", full.names = TRUE)
  # bio_rcp60 <- list.files(file.path(gsdms_data, "gcm_30s"), pattern = "*bc60*", full.names = TRUE)
  # bio_rcp85 <- list.files(file.path(gsdms_data, "gcm_30s"), pattern = "*bc85*", full.names = TRUE)



## >> Processing covariate layers ####
cov_files <- c(srtm, soil, popdens, landuse, bio_current, gcms)
all(lapply(cov_files, file.exists))

## step one_clip by e
infile <- cov_files
# file_in <- c(bioclim, srtm, soil, landuse)
outfile <- file.path(data_processed, paste0(tools::file_path_sans_ext(basename(infile)), "_cropped.", tools::file_ext(infile)))
reso_wgs <- res(readRDS(file.path(rdata_path, paste0("mask_", region, ".tif"))))
e <- extent(readRDS(file.path(rdata_path, paste0("mask_", region, ".tif"))))
mapply(gdalCrop, inpath = infile, outpath = outfile, MoreArgs = list(extent=e, resolution=reso_wgs, return = FALSE)) 

## step two_reproject: all to WGS84
infile <- outfile #list.files(data_processed, pattern = "_cropped", full.names = TRUE)
outfile <- sub("_cropped", "_aligned", infile)

new_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
new_res <- c(0.008333333,0.008333333)

  # if (region == "aus") {
  #   ## aus: Australian Albers
  #   new_crs = "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  #   new_res = c(...)
  # } else {
  #   ## vn: WGS84
  #   new_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  #   new_res = c(0.008333333,0.008333333)
  # }

file.res <-  ...

for (i in length(infile)) {
  gdalwarp(srcfile = infile[i], dstfile = outfile[i], s_srs = fileres[i], t_srs = new_crs, tr = new_res, verbose=TRUE)
}

mapply(gdalwarp, srcfile = infile, dstfile = outfile, MoreArgs = list(s_srs = wgs_crs, t_srs = new_crs, tr = reso, verbose=TRUE))

## step three_mask
## https://gitlab.unimelb.edu.au/garberj/gdalutilsaddons/-/blob/master/gdal_calc.R
outfile2 <- sub("_aligned", "_masked", outfile)
outfile2 <- paste0(tools::file_path_sans_ext(outfile2), ".tif")

mask_file <- file.path(rdata_path, paste0("mask_", region, ".tif"))
mapply(gdalmask, infile = outfile, mask = mask_file, outfile = outfile2, MoreArgs = list(output_Raster = FALSE, overwrite=FALSE, verbose=TRUE))
# file.remove(c(infile, outfile))


## Find min non-NA set values across mask and covariates and sync NAs ####
infile <- list.files(data_processed, pattern = "_masked", full.names = TRUE)
mask_file <- readRDS(file.path(rdata_path, paste0("mask_", region, ".tif")))
mask_file2 <- file.path(data_processed, paste0("mask_", region, "_nona", ".tif"))
file.copy(mask_file, mask_file2)

for(j in 1:length(infile)){
  print(paste0("processing file = ", j, " of ", length(infile) ,":", infile[j], " for updating mask"))
  # gdalmask(infile = mask_file, mask = infile[j], outfile = mask_file, output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
  gdalcalc(calc="((A>-9999) & (B==1))", infile = list(A=infile[j], B=mask_file2),outfile = mask_file2,
           NoDataValue=-9999, overwrite=TRUE)
}

for(j in 1:length(infile)){
  print(paste0("processing file = ", j, " of ", length(infile) ,":", infile[j], " using updated mask"))
  gdalmask(infile = infile[j], mask = mask_file2, outfile = infile[j], output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
}

## Check
summary(raster(mask_file))
summary(raster(mask_file2))
## jgarber comment: Im not sure what the summary file is doing, doesn't seem to count the NA's right, but if you subtract the two rasters and plot them, you will see that some of the cells in the continent have a difference of 1, meaning  that they were masked, but the no data values in the data files are now 0's in the new mask file, this should still mask them out when you run gdalmask, as it masks out every cell where the mask does not equal 1. So i think we are good to go!(see the photo below. White is 0, black is a nodatavalue)
freq(raster(mask_file))
freq(raster(mask_file2))
## jgarber comment: Ran the frequency and as you can see, it has removed some ones from the original mask and made them zeros, this should still work in the gdalmask function. From looking at the map, it seems like it is mostly picking up on large lakes? Looks like the great lakes, Lake victoria, and Lake Baikal are some of the removed cells.
gdalinfo(mask_file, stats = TRUE)
gdalinfo(mask_file2, stats = TRUE)

lapply(infile, gdalinfo, stats = TRUE)


# ## Extract cell-wise quartiles across GCM ####
# quartiles <- c("q1", "q2", "q3")
# 
# ## Files by models
# models <- tolower(models)
# mod_files <- paste0("gcm_", models)
# 
# for (i in 1:length(models)){
#   assign(mod_files[i], list.files(file.path(gsdms_data, "gcm_30s"), 
#                                   pattern = models[i], full.names = TRUE))
# }
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
# file.copy(list.files(gcm_quant_path, full.names = TRUE), file.path(rdata_path))
# rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
# gc()



## Bioregions layer (for Australia only) ####
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
writeRaster(bioreg_rast, file.path(rdata_path, "bioregions_aus.tif"))



## Sync NAs again?? ####

## >> Prepare GTAP data ####
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


## >> Biodiversity data ####
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


