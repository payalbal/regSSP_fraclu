rm(list = ls())
library("raster")
source(file.path(".", "R", "1_functions.R"))

countries <- c("vnm", "aus")
rcps <- c("1oC", "4oC")
scens <- c("26", "85")
j <- 1
timesteps <- c(1, 13, 23, 33, 43, 53) #these time steps we want outputs for to use in land use model
for (j in 1:length(countries)){
  country <- countries[j]
  landuse <- readRDS(file.path(".", "RData", paste0("covariates_", country, ".rds")))
  landuse <- landuse[[which(names(landuse) == "landuse")]]
  demand_2005 <- table(landuse[])[-c(7,8)]
  
  for (i in 1:length(rcps)){
    rcp <- rcps[i]
    yrs <- length(2018:2070)
    demand_2050 <- matrix(NA, nrow = yrs, ncol = length(demand_2005))
    demand_2050[1,] <- demand_2005
    
    #i) Urban land
    if (country == 'vnm'){
      urb_change <- (55739 - 33121.4)/33121.4 #Figures of urban population change from FAO
    } else if (country == 'aus') {
      urb_change <- (31346.1 - 21996.1)/21996.1 #Figures of urban population change from FAO
    }
    an_urb_change <- urb_change/yrs
    demand_2050[,1] <- round(demand_2005[1] * (1 + c(1:(yrs)) * c(0, rep(an_urb_change, yrs-1))))
    
    #ii) Cropland
    #Load proportions harvested in 2016
    total_bysector <- readRDS(file.path(".", "RData", paste0("harvested2016_", country, ".rds")))
    
    #Load projected sector landendowments
    m_agr <- readRDS(file = file.path(".","RData", paste0("gtap_landendowments_", country, "_", scens[i], ".rds")))
    
    #calculate mean sector output, weighted by fao estimates of harvests
    demand_2050[,2] <- round(demand_2005[2] * (1 + colMeans(total_bysector[,2] * m_agr[-which(rownames(m_agr)%in%c("frs", "ctl")),]/100 * length(total_bysector[,2]))))
    
    #wet and barren
    #demand_2050[,7] <- round(demand_2005[7])
    #demand_2050[,8] <- round(demand_2005[8])
    
    #Grass/shrub and forest demand are calculated during land-use simulations, based on their repsective mean suitabilities in the landscape, see 3_landuse.R
    
    demand <- demand_2050[timesteps,]
    
    colnames(demand) <- NULL
    print(rowSums(demand))
    
    saveRDS(demand, file.path(".", "output", paste0("landdemand_", country, "_", scens[i], ".rds")))
  }
}

