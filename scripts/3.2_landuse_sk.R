## Land use model for fractinal data
## Ref: https://github.com/kapitzas/ch2_landusemodel

## Set working environment ####

# setwd("./regSSP_fraclu/")

rm(list = ls())
x <- c("sp", "raster", "extraDistr")
lapply(x, require, character.only = TRUE)
source(file.path(".", "scripts", "0_functions_sklu.R"))
source(file.path(".", "scripts", "0_functions.R"))

regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' # on server - "./data"
rdata_path <- file.path(regSSP_birds_data, "RData") # on server - "./RData"
output_path <- file.path(regSSP_birds_data, "output") # on server - "./output"
files <- list.files(rdata_path, full.names = TRUE)


## Specify region ####
region = 'vn'
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))


## Load covariates ####
covariates_all <- readRDS(files[grepl(paste0("covariates_", region), files)])
landuse_layers <- grepl("copernicus", names(covariates_all))
landuse_layers <- names(covariates_all)[landuse_layers]
covs_exclude <- c("pa", "popd", landuse_layers)
covariates <- covariates_all[[-which(grepl(paste(covs_exclude, collapse = "|"), 
                                           names(covariates_all)))]] 


## Observed data: land use ####
lu_stack <- covariates_all[[landuse_layers]]
lu <- matrix(data = NA, ncell(reg_mask), nlayers(lu_stack))
for(i in 1: nlayers(lu_stack)){
  lu[,i] <- getValues(lu_stack[[i]])
}
lu <- na.omit(lu)

## Change classes as per landuse data
lu_classes <- c("urban", "crop", "forest", "grass", "other")
colnames(lu) <- lu_classes

## Checks: 
## Is total # cells - # of NA cells = dim of lu_mat_sub?
ncell(lu_stack[[1]]) - sum(is.na(getValues(lu_stack[[1]]))) == nrow(lu)
## Do all rows add up to 1?
## Note some minute residual decimal places, therefore rounding to 1
rs <- round(rowSums(lu, na.rm = TRUE, dims = 1)); all(rs==1)


## Convert data to data.frame #### 
reg_mask0 <- reg_mask
reg_mask0[which(is.na(reg_mask0[]))] <- 0
rpts <- rasterToPoints(reg_mask0, spatial=TRUE)
data <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
data <- data[,-1]
data <- na.omit(cbind(data, as.matrix(covariates)))

## Checks:
sum(!is.na(reg_mask[]))
nrow(data)


## Calculate neighbourhood statistics (mean values from moving window) ####
weights <- list(matrix(1, 3, 3, byrow= TRUE)) # size of window
weights <- rep(weights, length.out = 5)
ln <- neighbourhood(lu, cols = c(1 :5), weights, reg_mask) # ref functions script to see what's happening here
data <- cbind(data, ln)


## Correlation analysis and subset uncorr preds ####
preds <- colnames(correlations(data)) # correlation analysis and removes predictors from correlated pairs
data <- data[,colnames(data)%in%preds]


## Suitability model ####
## Specify model formula
form <- paste(preds, collapse="+") # linear terms

## Build model
subs <- sample(1:nrow(lu), 100000)
suitmod <- suitmodel(form = form, lu = lu[subs,], data = data[subs,], resolution = 1000, model = FALSE, maxit = 1000)

>>>> 
  ## Extract demand change from observed time steps (by country)
ts <- c(2015)
demand_all <- demand(landuse = lu, ts = ts ,inds = NULL, k = 5)[,1:(5+1)] #demand function pulls demand change from observed time steps and interpolates missing years


## Allocation algorithm ####
## Parameters for allocation algorithm: these can and probably need to be adjusted to achieve convergence and improve realism
allocpars <- list(
  stepsi = 1, # Stepsize for iterator. smaller values, slower convergence, higher chance for solution
  max_dev = 1, # Maximumn allowed deviation from prescribed demand.
  max_change = c(1,1,1,1,1), #by how much can land use grow in a cell? 1 means a land use can grow to occupy entire cell in one time step.
  min_change = - c(0.01, 0.05, 1, 1, 1), # opposite of above, i.e. urban and crop can only decrease by tiny fractions
  no_change = c(0.8, 0, 0, 0.05, 0), # growth threshold: in this case, when i.e. only when urban neighbourhood is > 0.8, urban can be allocated in a cell that currently has 0 urban. Prevents i.e. urban growth in the middle of a forest.
  suit_adj = 1e-14, # Amount of landcape-wide adjustment of predicted suitability when demand allocations increase the difference between allocated and prescribed supply.
  resolution = 10000
)

ts <- c(2000, 2006, 2012, 2018)
lu_out <- lu
lu_ts <- list()

# 3.b Iterate through time steps
i <- 1
for(i in 1:(length(ts)-1)){
  cat(paste0('\n', "Predicting suitability for time step ", i))
  
  # 3.c Caluclate neighbourhood maps for current time step
  ln <- neighbourhood(lu_out, 1:5, weights, reg_mask)
  
  #3.d Predict the suitability model, adding the neigbourhood maps as covariates
  sm <- predict(suitmod, newdata = cbind(dat, ln), type = "probs")
  
  lu_pred <- matrix(NA, nrow(lu), ncol(lu))
  
  # 3.d get demand change for time step
  dmd <- demand_all #j instead of 12
  dmd <- dmd[which(dmd[,1]%in%ts),-1]
  dmd_t0 <- dmd[i,]
  
  if(i > 1){
    dmd_t0 <- colMeans(lu_ts[[i-1]])
  }
  dmd_ts <- rbind(dmd_t0, dmd[i+1,])
  
  # 3.e Allocate, see functions script for what this does
  lu_pred <- allocation(lu = lu_out,
                        ln = ln, 
                        sm = sm, 
                        params = allocpars, 
                        dmd = dmd_ts)
  
  
  cat('\n')
  # 3.f Store predicted maps
  lu_ts[[i]] <- lu_out <- lu_pred
}

test <- reg_mask
test[!is.na(test[])] <- lu_ts[[3]][,1]
plot(test)

