## SDMs

## Set up work environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

# setwd("./regtrade/")
x <- c('sp', 'raster', "data.table", 'doParallel', "dismo", "rJava")
lapply(x, require, character.only = TRUE)
options( java.parameters = "-Xmx48g" )
# some info if java doens't initialise correctly
#https://support.apple.com/kb/DL1572?viewlocale=en_US&locale=en_US
#Error in .jnew("RectangularArrayBuilder", .jcast(array), dim) : 
#  java.lang.OutOfMemoryError: Java heap space
#to load rJava start RStudio from shell, using: LD_LIBRARY_PATH=$(/usr/libexec/java_home)/jre/lib/server: open -a 
#options("java.home"="/Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/jre")
source(file.path(".","R", "1_functions.R"))


regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' 
## change as per server: boab = "./data"
rdata_path <- file.path(regSSP_birds_data, "RData") 
## change as per server: boab - "./RData"
output_path <- file.path(regSSP_birds_data, "output") 
## change as per server: boab - "./output"


## Define global parameters ####
regions <- c("til", "aus")
bias_preds <- c("dibu", "diro", "pa", "popd", "roughness")
not_used <- c("srtm", "popd", "dico", "carb", "nitro", "awco", "bulk", "landuse") #these are not used in SDM
files <- list.files(rdata_path, full.names = T)
# not_used <- c(not_used, "soilnitro") # might be used later? not in covariates.


## I. Data preparation ####
for(region in regions){
  
  ## Load data ####
  covariates <- files[grepl(paste0("covariates_", region), files)]
  covariates <- readRDS(covariates)
  
  ## Create bias layer to account for sampling bias ####
  ## Bias layer predicted based on bias_preds
  birds <- readRDS(file.path(rdata_path, paste0("occ_", region, ".rds")))
  birds <- as.data.frame(birds)
  colnames(birds)[1:2] <- c("long", "lat")
  subs <- sample(1:nrow(birds), size = 5000)
  obs <- birds[subs, c(1,2)]
  covariates_bias <- covariates[[which(names(covariates)%in%bias_preds)]]
  print(paste0("Modelling bias layer for region ", region, " ..."))
  mod <- dismo::maxent(x= covariates_bias, p= obs,nbg = 10000, factors = "pa", args= c("randomtestpoints=25","-J", "-P", "-p", "-h", "threshold=FALSE"))
  print(paste0("Predicting bias layer for region ", region, " ..."))
  bias <- dismo::predict(mod, covariates_bias)
  saveRDS(readAll(bias), file = file.path(output_path, paste0("bias_", region, ".rds")))
  rm(birds, subs, obs, covariates_bias, mod, bias)
  
  
  ## Discard correlated covariates ####
  covariates_subset <- covariates[[-which(names(covariates)%in%c(not_used, bias_preds))]]
  cov_df <- getValues(covariates_subset)
  cov_df <- na.omit(cov_df)
  preds <- colnames(correlations(cov_df))
  saveRDS(preds, file = file.path(output_path, paste0("preds_", region, ".rds")))
  rm(covariates_subset, cov_df, subs, preds)
  
  
  ## Australia-specific: Bioregion layer ####
  ## Create list of adjacent bioregions for each unique bioregion (for background sampling)
  if(region == "aus"){
    bio_regs <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
    ecoregs <- na.omit(unique(bio_regs[]))
    ecoregs <- sort(ecoregs)
    adjs_list <- list()
    for (i in ecoregs){
      adjs <- adjacent(bio_regs, which(bio_regs[] == i), pairs = F, include = T)
      adjs_list[[i]] <- unique(bio_regs[adjs], na.rm = T)
      print(i)
    }
    saveRDS(adjs_list, file = file.path(output_path, "adjs_list.rds"))
  }
  rm(bio_regs, ecoregs, adjs_list, i, adjs)
  
  
  ## Sample background points using bias layer ####
  if (region == "til"){
    bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))
    inds <- which(!is.na(bias_rast[]))
    probs <- round(bias_rast[inds],3)
    inds_all <- sample(inds, p = probs, size = 10000)
    saveRDS(inds_all, file = file.path(output_path, paste0("bgp_", region, "_bias.rds")))
  }
  rm(bias_rast, inds, probs, inds_all)
  
  ## For aus, background points are constrained to the species bioregion plus one adjacent
  if(region == "aus"){
    bio_regs <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
    bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))
    ## Min set non-NA values
    bias_rast <- mask(bias_rast, bio_regs)
    bio_regs <- mask(bio_regs, bias_rast)
    probs <- round(bias_rast[],3)
    aves_obs <- readRDS(file.path(rdata_path, paste0("occ_",region,".rds")))
    
    species <- sort(unique(aves_obs$species))
    biases <- c(T, F)
    adj_list <- readRDS(file.path(output_path, "adjs_list.rds"))
    for (bias in biases){
      samples_list <- list()
      for (i in 1:length(species)){
        subs <- aves_obs[as.character(aves_obs$species)%in%species[i],]
        spec_bio_regs <- na.omit(unique(extract(bio_regs, subs[,c(1,2)])))
        adj_bioregs <- na.omit(unique(unlist(adj_list[spec_bio_regs])))
        inds <- which(bio_regs[]%in%adj_bioregs)
        sample_size <- 10000
        if (length(inds) > 1000000){
          inds <- sample(inds, size = 1000000)
        }
        if (length(inds) < 10000){
          sample_size <- length(inds)
        }
        if(bias == T){
          p <- round(probs[inds], 3)
          samples <- sample(x=inds, size = sample_size, p = p, replace = T)
        }else if(bias == F){
          samples <- sample(x=inds, size = 10000)
        }
        samples_list[[i]] <- samples
        print(i)
      }
      if(bias == T){
        saveRDS(samples_list, file = file.path(output_path, "bgp_aus_bias_constrained.rds"))
      }else if(bias == F){
        saveRDS(samples_list, file = file.path(output_path, "bgp_aus_random_constrained.rds"))
      }
    }
  }
}

rm(list=ls())
gc()



## II. Models ####
## ***** restart here ***** ####

## Set work environment ####
# setwd("./regtrade/")
x <- c('sp', 'raster', "data.table", 'doParallel', "dismo", "rJava", 'foreach', 'iterators')
lapply(x, require, character.only = TRUE)
options( java.parameters = "-Xmx48g" )
# some info if java doens't initialise correctly
#https://support.apple.com/kb/DL1572?viewlocale=en_US&locale=en_US
#Error in .jnew("RectangularArrayBuilder", .jcast(array), dim) : 
#  java.lang.OutOfMemoryError: Java heap space
#to load rJava start RStudio from shell, using: LD_LIBRARY_PATH=$(/usr/libexec/java_home)/jre/lib/server: open -a 
#options("java.home"="/Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/jre")
source(file.path(".","R", "1_functions.R"))


## Maxent trobleshooting
maxent()
#This is only if maxent() doens't work because java could not be found.
# system("ps")
# system("pkill -f R")
# #If Maxent doesn't work, do this:
# remove.packages("dismo")
# remove.packages("rJava")
# quit(save = "no", status = 0, runLast = F)
# install.packages("dismo")
# install.packages("rJava")
# quit(save = "no", status = 0, runLast = TRUE)
# file.copy("~/maxent.jar", "~/.r-dir/R/library/dismo/java/maxent.jar") #maxent.3.3.3
# file.exists("~/.r-dir/R/library/dismo/java/maxent.jar")
maxent()


regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' 
## change as per server: boab = "./data"
rdata_path <- file.path(regSSP_birds_data, "RData") 
## change as per server: boab - "./RData"
output_path <- file.path(regSSP_birds_data, "output") 
## change as per server: boab - "./output"

## Global parameters ####
regions <- c("til", "aus")
ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))
treatments <- c("pre", "ind", "dir", "agg")

## Aus results only...
# scens <- scens[c(1,3,7)]
# scens_rcps <- scens_rcps[c(1,3,7)]

## ***** Specify region for analysis ***** ####
##  Rerun for c("til", "aus")
region <- 'til'


## Load data ####
## Observation data
obs <- readRDS(file.path(rdata_path, paste0("occ_",region,".rds")))
species <- sort(as.character(unique(obs$species)))

## Covariates for prediction, write into temp folder structure
##  Note: Covariates used for model building and predictions 
##    are different for Vietnam, see Methods
preds <- readRDS(file.path(output_path, paste0("preds_",region,".rds")))

if(region == "til"){
  region_preds <- "vn"
}else{
  region_preds <- region
}

temp <- paste0(file.path(output_path, "temp"))
type <- c("covariates", "predlu", "bio")
unlink(temp, recursive = T)
for(k in 1:3){
  if(type[k]%in%type[c(2,3)]){
    for(i in 1:length(scens)){
      if(type[k] == "predlu"){
        path <- output_path
        pattern <- paste0(type[k], "_", scens[i], "_", region_preds)
      }
      if(type[k] == "bio"){
        path <- rdata_path
        pattern <- paste0(type[k], scens_rcps[i],"_", region_preds)
      }
      covariates <- readRDS(list.files(path, pattern = pattern, full.names = T))
      covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
      writeToDisk(covariates, folder = file.path(temp, pattern))
    }
  }else{
    covariates <- readRDS(list.files(rdata_path, pattern = paste0("covariates_", region_preds), full.names = T))
    covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
    writeToDisk(covariates, folder = file.path(temp, "pres"))
  }
}

## Covariates for modelling, write into temp folder structure
covariates <- readRDS(list.files(rdata_path, pattern = paste0("covariates_", region), full.names = T))
covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
writeToDisk(covariates, folder = file.path(temp, "models"))
covariates <- stack(list.files(file.path(temp, "models"), full.names = T))

## Bias layer and list of predictors to retain
bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))

## Background points
if(region == "aus"){
  bio_region <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
  samples_list <- readRDS(file.path(output_path, "bgp_aus_bias_constrained.rds"))
} else {
  inds <- readRDS(file.path(output_path, paste0("bgp_", region, "_bias.rds")))
  bgp <- SpatialPoints(xyFromCell(bias_rast, cell = inds))
  bg <- extract(covariates, bgp)
}

## Specify model parameters
output <- file.path(output_path, paste0("models_", region))
if(!dir.exists(output)){dir.create(output)}
factors <- "landuse"
logfile <- file.path(output, paste0(region, "_log.txt"))
writeLines(c(""), logfile)
l <- length(species)
predictions <- T


## Build Maxent model ####
## Load Cluster
cl <- makeCluster(12)
registerDoParallel(cl)
cl <- makeCluster(12)
registerDoParallel(cl)
clusterCall(cl, function() library(sp))
clusterCall(cl, function() library(raster))
clusterCall(cl, function() library(data.table))
clusterCall(cl, function() library(dismo))
clusterCall(cl, function() library(rJava))
clusterCall(cl, function() library(foreach))
clusterCall(cl, function() library(iterators))
clusterCall(cl, function() library(parallel))
clusterCall(cl, function() library(doParallel))


results <- foreach(i = 1:l, 
                   .packages = c('sp', 'raster', "data.table", 'doParallel', 
                                 "dismo", "rJava",'foreach', 'iterators')) %dopar% {
  cov <- covariates
  cat(paste("Starting model",i,"\n"), file = logfile, append = T)
  subs <- obs[obs$species == species[i],]
  
  ## Load background data
  if(region == "aus"){
    inds <- samples_list[[i]]
    bgp <- SpatialPoints(xyFromCell(bias_rast, cell = inds))
    bg <- extract(cov, bgp)
  }
  
  ## Prepare data for model building
  sp <- SpatialPoints(subs[,c(1,2)])
  pr <- extract(cov, sp)
  occ <- c(rep(0, 1, nrow(bg)), rep(1, 1, nrow(pr)))
  
  cov_df <- data.frame(rbind(bg, pr))
  
  ## Determine covariate importance
  ## Initial model run
  mod <- tryCatch(maxent(cov_df, 
                         p= occ, 
                         factors = factors, 
                         args= c("randomtestpoints=25",
                                 "-J", "-p", "-h", 
                                 "threshold=FALSE")),
                  error = function(e) NA)
  
  ## Determine covariates with perm importance < 1 and remove from covariate set (see Methods)
  if(class(mod) != "MaxEnt"){
    return(NA)
  } else {
    
    res <- data.frame("names" = as.character(rownames(mod@results)), "results" = as.numeric(mod@results))
    perm <- res[grep("permutation.importance", res$names),]
    kickout <- perm$names[which(perm$results < 1)]
    kickout <- gsub('.permutation.importance', '', kickout)
    
    if(length(kickout) > 0){
      cov_df <- cov_df[,-which(names(cov)%in%kickout)]
      names.new <- colnames(cov_df)
    }else{
      cov_df <- cov_df
      names.new <- colnames(cov_df)
    }
    
    factors.new <- factors[factors%in%names(cov_df)]
    msg1 <- paste0("landuse covariate discarded for species i = ", i, " ", species[i], 
           " for region = ", region, ": ", length(factors.new) == 0)
    msg1
    cat(msg1, file = logfile, append = T)
    
    ## Cross-validated model
    mod.new <- tryCatch(dismo::maxent(cov_df, 
                                      p= occ,
                                      factors = factors.new, 
                                      args= c("replicatetype=crossvalidate", 
                                              "replicates=4", 
                                              "-J","-P", "-p", "-h", 
                                              "threshold=FALSE")), 
                        error = function(e) NA)
    
    modresults <- mod.new@results
    rm(mod, mod.new)
    gc()
    
    ## Final model
    out_spec <- file.path(output, species[i])
    mod.final <- tryCatch(dismo::maxent(cov_df, 
                                        p= occ,
                                        factors = factors.new,
                                        args= c("randomtestpoints=25",
                                                "-J","-P", "-p", "-h", 
                                                "threshold=FALSE", "writeplotdata=TRUE"), path = out_spec), 
                          error = function(e) NA)
    
    
    if(class(mod.final) != "MaxEnt"){
      return(NA)
    } else {
      if(predictions == F) { next }
      

      ## Prediction ####
      cat(paste("Starting prediction",i,"\n"), file = logfile, append = T)
      
      ##Get maxSSS threshold to deliniate suitable habitat (Liu et al 2016) 
      thresh_maxSSS <- mod.final@results[names(mod.final@results[,1])%in%"Maximum.test.sensitivity.plus.specificity.Cloglog.threshold"]
      ## Note: Previously 'Maximum.test.sensitivity.plus.specificity.logistic.threshold' which is not listed anymore...???
      
      ## Get present cov set and static variables (don't change) for predictions
      cov <- stack(list.files(file.path(temp, "pres"), full.names = T))
      sta <- cov[[-which(grepl(paste0(c("landuse", "bio"), collapse = "|"), names(cov)))]]
      area_maxSSS <- numeric()
      cells_maxSSS <- list()
      m <- 0 # index counter for cells_maxSSS
      
      for (k in 1:length(treatments)) {
        
        for (j in 1:length(scens)) { #scens are a combination of RCP * Quartile
          m <- m+1
          ##Load data to predict under different treatments and scenarios
            ## Present prediction (current land use and climate)
          if (treatments[k] == "pre") {
            dyn <- cov[[which(grepl("bio", names(cov)))]]
            luc <- cov[[which(grepl("landuse", names(cov)))]]
            
            ##'Aggregated' Land use predcitors as per ssp and climate predictors as per rcps
          } else if(treatments[k] == "agg") {
            dyn <- stack(list.files(file.path(temp, paste0("bio", scens_rcps[j], "_", region_preds)), full.names = T))
            luc <- stack(list.files(file.path(temp, paste0("predlu_", scens[j], "_", region_preds)), full.names = T))
            
            ## Direct: Land use stays the same (current), only clime predictors change as per rcps
          } else if(treatments[k] == "dir") {
            dyn <- stack(list.files(file.path(temp, paste0("bio", scens_rcps[j], "_", region_preds)), full.names = T))
            luc <- cov[[which(grepl("landuse", names(cov)))]]

            ## Indirect: Land use changes as per ssp, rcp predictors stay the same (current)
          } else if(treatments[k] == "ind") {
            dyn <- cov[[which(grepl("bio", names(cov)))]]
            luc <- stack(list.files(file.path(temp, paste0("predlu_", scens[j], "_", region_preds)), full.names = T))
          }
          
          ## Combine static, dynamic and landuse covariates
          cov.new <- stack(sta, dyn, luc)
          cov.new <- stack(cov.new[[which(names(cov.new)%in%names.new)]])
          
          ## Predicted map as per region
          map <- dismo::predict(mod.final, cov.new)
          rm(cov.new)
          
          ## Apply bioregion constraint for aus (see Methods)
          if (regions == "aus"){
            map[which(!bio_region[]%in%na.omit(unique(bio_region[inds])))] <- NA
          }
          
          ## Get number of cells above or equal to maxSSS threshold: This is deemed suitable habitat
          area_maxSSS <- c(area_maxSSS, length(which(map[] >= thresh_maxSSS)))
          if (grepl("q2", scens[i])){
          cells_maxSSS[[m]] <- which(map[] >= thresh_maxSSS)}
          if(treatments[k] == "pre") { break }
        }
      }
    }
  }
  
  scens
  rm(dyn, luc, cov.new, mod.final)
  gc()
  list(modresults, area_maxSSS, i, species[i], cells_maxSSS)
}
saveRDS(results, file = file.path(output_path, paste0("results_", region, ".rds")))
unlink(temp, recursive = T)
stopCluster(cl)

## Structure of results
## results is a list of length(species)
## results[[x]][[1]] = maxent model output
## results[[x]][[2]] = areas_maxSSS vector of length = 1 pres treatment + length(treatments) -1 * length(scens)
## results[[x]][[3]] = index of species
## results[[x]][[4]] = name of species
## results[[x]][[5]] = list of cells_maxSSS vectors of length = 1 pres treatment + length(treatments) -1 * length(scens)

## ***** Rerun scipt up till here for all regions ***** ####


## Subset results for plotting
res <- results
temp <- list()

for (i in 1:length(res)){
  temp2 <- list(res[[i]][[5]][[1]], res[[i]][[5]][[23]], res[[i]][[5]][[24]], res[[i]][[5]][[25]])
  temp[[i]] <- list(res[[i]][[1]], res[[i]][[2]], res[[i]][[3]], res[[i]][[4]], temp2)
}

saveRDS(temp, file = file.path(output_path, paste0("results2_", region, ".rds")))