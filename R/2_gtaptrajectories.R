rm(list = ls())
library(sp)
library(raster)
source(file.path(".", "R", "1_functions.R"))


rcps <- c("26", "85")
timesteps <- c(1, 13, 23, 33, 43, 53) # these time steps we want outputs for to use in land use model
lu_names <- c("urban", "crop", "herb", "shrub", "Oforest", "Cforest")
year_range <- 2018:2070

for (region in c("vnm", "aus")){
  landuse <- readRDS(file.path(".", "RData", paste0("covariates_", region, ".rds")))
  landuse <- landuse[[which(names(landuse) == "landuse")]]
  lu_yr0 <- table(landuse[])[-c(7,8)] # lu for base year
  
  for (j in 1:length(rcps)){
    yrs <- length(year_range)
    lu_byYear <- matrix(NA, nrow = yrs, ncol = length(lu_yr0)) #empty matrix to store lu change trajaectories
    colnames(lu_byYear) <- lu_names
    rownames(lu_byYear) <- paste0("yr", year_range)
    lu_byYear[1,] <- lu_yr0
    
    ## Urban
    ## Data source: http://www.fao.org/faostat/en/#data/LC
    urban_cover <- as.data.frame(read.csv("./raw_data/fao_data/fao_urban_landcover.csv", header = T))
    urban_cover <- urban_cover[,c("Area", "Year", "Value")]
    urban_cover$Area <- gsub("Australia", "aus", urban_cover$Area)
    urban_cover$Area <- gsub("Viet Nam", "vnm", urban_cover$Area)
    
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
    #Load proportions harvested in 2016
    total_bysector <- readRDS(file.path(".", "RData", paste0("harvested2016_", region, ".rds")))
    
    #Load projected sector landendowments
    m_agr <- readRDS(file = file.path(".","RData", paste0("gtap_landendowments_", region, "_", rcps[j], ".rds")))
    
    #calculate mean sector output, weighted by fao estimates of harvests
    lu_byYear[,2] <- round(lu_yr0[2] * (1 + colMeans(total_bysector[,2] * m_agr[-which(rownames(m_agr)%in%c("frs", "ctl")),]/100 * length(total_bysector[,2]))))
    
    ## Rest of the lu classes (herb, shrub, Oforest, Cforest) are calculated during land-use simulations, based on their repsective mean suitabilities in the landscape, see 3_landuse.R
    
    demand <- lu_byYear[timesteps,]
    
    colnames(demand) <- NULL
    print(rowSums(demand))
    
    saveRDS(demand, file.path(".", "lu_output", paste0("landdemand_", region, "_", rcps[j], ".rds")))
  }
}




urb_change <- (55739 - 33121.4)/33121.4
} else if (region == 'aus') {
  urb_change <- (31346.1 - 21996.1)/21996.1
}
an_urb_change <- urb_change/yrs
