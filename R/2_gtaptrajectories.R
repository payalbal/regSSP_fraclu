## Estimate land use demands for urban and crop lu classes using FAO data
##
## Land endowments for other lu classes (herb, shrub, Oforest, Cforest) are estimated based on land-use simulations (accounting for mean suitability of land for LU class) see 3_landuse.R


## Set up work environment
rm(list = ls())
library(sp)
library(raster)
source(file.path(".", "R", "1_functions.R"))
regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' # change as per server: boab = "./data"
layer_path <- file.path(regSSP_birds_data, "RData") # change as per server: boab - "./RDdata"
out_path <- file.path(regSSP_birds_data, "lu_output")


## Define global variables
ssps <- paste0("ssp", 1:3)
timesteps <- c(1, 11, 21, 31, 41, 51) # these time steps we want outputs for to use in land use model
lu_names <- c("urban", "crop", "herb", "shrub", "Oforest", "Cforest")
year_range <- 2019:2070


## Estimate land demand for lu classes for each timestep
for (region in c("vn", "aus")){
  print(paste0("processing ", region, "... "))
  landuse <- readRDS(file.path(layer_path, paste0("covariates_", region, ".rds")))
  landuse <- landuse[[which(names(landuse) == "landuse")]]
  lu_yr0 <- table(landuse[])[-c(7,8)] # lu for base year (excluding classes 7 & 8)
  names(lu_yr0) <- lu_names
  
  for (j in 1:length(ssps)){
    yrs <- length(year_range)
    lu_byYear <- matrix(NA, nrow = yrs, ncol = length(lu_yr0)) #empty matrix to store lu change trajaectories
    colnames(lu_byYear) <- lu_names
    rownames(lu_byYear) <- paste0("yr", year_range)
    lu_byYear[1,] <- lu_yr0
    
    ## Urban
    ## Data source: http://www.fao.org/faostat/en/#data/LC
    urban_cover <- as.data.frame(read.csv(file.path(regSSP_birds_data, "fao_data", "fao_urban_landcover.csv"), header = T))
    urban_cover <- urban_cover[,c("Area", "Year", "Value")]
    urban_cover$Area <- gsub("Australia", "aus", urban_cover$Area)
    urban_cover$Area <- gsub("Viet Nam", "vn", urban_cover$Area)
    
    reg_dat <- urban_cover[which(urban_cover$Area == region),]
    nT <- reg_dat[which(reg_dat$Year == max(reg_dat$Year)),]$Value
    n0 <- reg_dat[which(reg_dat$Year == min(reg_dat$Year)),]$Value
    t <- nrow(reg_dat)
    
    ## Estimate rate of proportional change using discrete growth model
    ## nt = r^t * n0
    ## therefore, r = (nt/n0)^1/t
    urb_ratePropChange <- (nT/n0)^(1/t)
    ## or exp((log(nT) - log(n0)) * (1 / t))
    
    ## OR Estimate asoulte rate of change scaled by change in year T
    # an_urb_change <- ((nT-n0)/nT)*(1/t)
    
    ## Estimate area for all years, using discrete growth model
    for(i in 2:yrs){
      lu_byYear[i,1] <- lu_byYear[i-1,1] * urb_ratePropChange
    }
    lu_byYear[,1] <- round(lu_byYear[,1])
    ## OR if using an_urb_change:
    # lu_byYear[,1] <- round(lu_yr0[1] * (1 + c(1:(yrs)) * c(0, rep(an_urb_change, yrs-1))))
    
    ## Crop
    ## Load proportions harvested in yr0 (here, as 2019)
    yr0_bysector <- readRDS(file.path(layer_path, paste0("harvested2016_", region, ".rds")))
    
    ## Load estimated land endowments by sector for 2020 - 2070 (estimated from gtap data)
    yrsAll_bysector <- readRDS(file = file.path(layer_path, paste0("gtap_landendowments_", region, "_", ssps[j], ".rds")))
    
    ## Calculate mean sector output, weighted by fao estimates of harvests
    lu_byYear[,2] <- round(lu_yr0[2] * (1 + colMeans(yr0_bysector[,2] * yrsAll_bysector[which(rownames(yrsAll_bysector)%in% yr0_bysector$gtap_sector),]/100 * nrow(yr0_bysector))))
    
    demand <- lu_byYear[timesteps,]
    colnames(demand) <- NULL
    print(paste0("Land demand for ", region, " & ", ssps, ":")[j])
    print(rowSums(demand))
    print("----------------------")
    saveRDS(demand, file.path(out_path, paste0("landdemand_", region, "_", ssps[j], ".rds")))
  }
}




# urb_change <- (55739 - 33121.4)/33121.4
# } else if (region == 'aus') {
#   urb_change <- (31346.1 - 21996.1)/21996.1
# }
# an_urb_change <- urb_change/yrs
