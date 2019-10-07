rm(list = ls())
#required packages
require("raster")
library("data.table")
library("raster")
library("doParallel")
library("dismo")
library("rJava")
library("sp")

source(file.path(".","R", "1_functions.R"))

#-------------------------------#
#####---I. DATA PREP FOR SDMS####
#-------------------------------#

#1.a Set some global parameters
country_abbr <- "aus" #aus or til (til is tile29 of the worldclim data, see Methods)
layer_path <- file.path(".", "RData")
output_path <- file.path(".", "output")
bias_preds <- c("dibu", "diro", "pa", "popd", "roughness")
not_used <- c("srtm", "popd", "dico", "carb", "nitro", "awco", "bulk", "soilnitro", "landuse") #these are not used in SDM

options( java.parameters = "-Xmx48g" )

#1.b Make bias layers

# some info if java doens't initialise correctly
#https://support.apple.com/kb/DL1572?viewlocale=en_US&locale=en_US
#Error in .jnew("RectangularArrayBuilder", .jcast(array), dim) : 
#  java.lang.OutOfMemoryError: Java heap space
#to load rJava start RStudio from shell, using: LD_LIBRARY_PATH=$(/usr/libexec/java_home)/jre/lib/server: open -a 
#options("java.home"="/Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/jre")

countries <- c("aus", "til")
for(i in 1:2){
  country_abbr <- countries[i]
  files <- list.files(file.path(".", "RData"), pattern = country_abbr, full.names = T)
  files <- readRDS(files[grepl(paste0("covariates_", country_abbr), files)])
  files <- files[[which(names(files)%in%bias_preds)]]
  birds <- readRDS(file.path(".", "RData", paste0("occ_", country_abbr, ".rds")))
  subs <- sample(1:nrow(birds), size = 5000)
  obs <- birds[subs ,c(1,2)]
  
  #Only need spp names & lat/long cols
  print(paste0("Modelling", country_abbr))
  mod <- dismo::maxent(x= files, p= obs,nbg = 10000, factors = "pa", args= c("randomtestpoints=25","-J", "-P", "-p", "-h", "threshold=FALSE"))
  print(paste0("Predicting", country_abbr))
  bias <- dismo::predict(mod, files)
  saveRDS(readAll(bias), file = file.path(output_path, paste0("bias_", country_abbr, ".rds")))
}

#1.b) Correllation anlaysis
files <- list.files(layer_path, full.names = T)
covariates <- readRDS(files[which(grepl(paste0("covariates_", country_abbr), files))])
covariates <- covariates[[-which(names(covariates)%in%c(not_used, bias_preds))]]
cov_df <- getValues(covariates)
cov_df <- na.omit(cov_df)
subs <- sample(1:nrow(cov_df), size = 10000)

preds <- rownames(reduce_predset(cor(cov_df[subs,])))
saveRDS(preds, file = file.path(output_path, paste0("preds_", country_abbr, ".rds")))

#1.c) Australia-specific: determine adjacent bioregions (for background sampling)
if(country_abbr == "aus"){
  bio_regs <- readRDS(file.path(".", "RData", "bioregions_aus.rds"))
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

#1.d) Sampling background points using bias layer
if (country_abbr == "til"){
  dens_rast <- readRDS(file.path(output_path, paste0("bias_", country_abbr, ".rds")))
  inds <- which(!is.na(dens_rast[]))
  probs <- round(dens_rast[inds],3)
  inds_all <- sample(inds, p = probs, size = 10000)
  saveRDS(inds_all, file = file.path(".", "output", paste0("bgp_", country_abbr, "_bias.rds")))
}

if(country_abbr == "aus"){
  #- these backgorund points are constrained to the species bioregion plus one adjacent
  country_abbr <- "aus" #only australia
  bio_regs <- readRDS(file.path(layer_path, "bioregions_aus.rds"))
  dens_rast <- readRDS(file.path(output_path, paste0("bias_", country_abbr, ".rds")))
  dens_rast <- mask(dens_rast, bio_regs)
  bio_regs <- mask(bio_regs, dens_rast)
  probs <- round(dens_rast[],3)
  aves_obs <- readRDS(file.path(layer_path, paste0("occ_",country_abbr,".rds")))
  
  species <- sort(unique(aves_obs$species))
  biases <- c(T, F)
  adj_list <- readRDS(file.path(output_path, "adjs_list.rds"))
  for (j in 1:2){
    samples_list <- list()
    bias <- biases[j]
    for (i in 1:length(species)){
      subs <- aves_obs[as.character(aves_obs$species)%in%species[i],]
      spec_bio_regs <- na.omit(unique(extract(bio_regs, subs[,c(1,2)])))
      adj_bioregs <- na.omit(unique(unlist(adjs_list[spec_bio_regs])))
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

#-------------------#
#2. LOAD SDM DATA####
#-------------------#
#2.a Global parameters
country_abbr <- "til" #aus or til, for australia and bioclim tile29, which encompasses Vientam
scens <- c(paste0("q1", c("_26", "_85")), paste0("q2", c("_26", "_85")), paste0("q3", c("_26", "_85")))
treatments <- c("pre", "ind", "dir", "agg")

#2.b) Bias layer and predictors we wish to retain
dens_rast <- readRDS(file.path(layer_path, paste0("bias_", country_abbr, ".rds")))
preds <- readRDS(file.path(output_path, paste0("preds_",country_abbr,".rds")))

#2.c) Observation data
obs <- readRDS(file.path(layer_path, paste0("occ_",country_abbr,".rds")))
species <- sort(as.character(unique(obs$species)))

#2.d) Covariates

#i) assign country id for predictions - model building and predictions are different for Vietnam, see Methods
if(country_abbr == "til"){
  country_abbr_preds <- "vnm"
}else{
  country_abbr_preds <- country_abbr
}

#ii) load layers from rds files and write into folder structure, so they are stored locally and we don't have to load them into memory
temp <- paste0(file.path(output_path, "temp"))
type <- c("covariates", "landuse", "bio")
unlink(temp, recursive = T)
k <- i <- 2
for(k in 1:3){
  if(type[k]%in%type[c(2,3)]){
    for(i in 1:length(scens)){
      pattern <- paste0(type[k], scens[i],"_", country_abbr_preds)
      if(type[k] == "landuse"){
        path <- output_path
      }
      if(type[k] == "bio"){
        path <- layer_path
      }
      covariates <- readRDS(list.files(path, pattern = pattern, full.names = T))
      covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
      writeToDisk(covariates, folder = file.path(temp, pattern))
    }
  }else{
    covariates <- readRDS(list.files(layer_path, pattern = paste0("covariates_", country_abbr_preds), full.names = T))
    covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
    writeToDisk(covariates, folder = file.path(temp, "pres"))
  }
}

#iii) Load covariates for modelling, write into temp folder structure
covariates <- readRDS(list.files(layer_path, pattern = paste0("covariates_", country_abbr), full.names = T))
covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
writeToDisk(covariates, folder = file.path(temp, "models"))
covariates <- stack(list.files(file.path(temp, "models"), full.names = T))

#2.d) Get Background points from file
if(country_abbr == "aus"){
  bio_region <- readRDS(file.path(layer_path, "bioregions_aus.rds"))
  samples_list <- readRDS(file.path(output_path, "bgp_aus_bias_constrained.rds"))
} else {
  inds <- readRDS(file.path(output_path, paste0("bgp_", country_abbr, "_bias.rds")))
  bgp <- SpatialPoints(xyFromCell(dens_rast, cell = inds))
  bg <- extract(covariates, bgp)
}

#2.e) Model parameters
output <- file.path(output_path, paste0("models_", country_abbr))
factors <- "landuse"
logfile <- file.path(output, paste0(country_abbr, "_log.txt"))
writeLines(c(""), logfile)
l <- length(species)
predictions <- T

#---------------#
#3. BUILD SDM####
#---------------#

#3.a) Load Cluster
cl <- makeCluster(12)
registerDoParallel(cl)

results <- foreach(i = 1:l, .packages = c("dismo")) %dopar% {
  cov <- covariates
  cat(paste("Starting model",i,"\n"), file = logfile, append = T)
  subs <- obs[obs$species == species[i],]
  
  #3.b) Load background data
  if(country_abbr == "aus"){
    inds <- samples_list[[i]]
    bgp <- SpatialPoints(xyFromCell(dens_rast, cell = inds))
    bg <- extract(cov, bgp)
  }
  
  #3.c) Prep data for model building
  sp <- SpatialPoints(subs[,c(1,2)])
  pr <- extract(cov, sp)
  occ <- c(rep(0, 1, nrow(bg)), rep(1, 1, nrow(pr)))
  
  cov_df <- data.frame(rbind(bg, pr))
  
  #3.d) Determine covariate importance
  
  #i) Initial model run
  mod <- tryCatch(maxent(cov_df, 
                         p= occ, 
                         factors = factors, 
                         args= c("randomtestpoints=25",
                                 "-J", "-p", "-h", 
                                 "threshold=FALSE")),
                  error = function(e) NA)
  
  #ii) Determine covariates with perm importance < 1 and remove from covariate set (see Methods)
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
    
    #3.e) Cross-validated model
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
    
    #3. f) Final model
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
      
      #-----------------#
      #4. PREDICTIONS####
      #-----------------#
      
      cat(paste("Starting prediction",i,"\n"), file = logfile, append = T)
      
      #4. a Get maxSSS threshold to deliniate suitable habitat 
      thresh_maxSSS <- mod.final@results[names(mod.final@results[,1])%in%"Maximum.test.sensitivity.plus.specificity.logistic.threshold"]
      
      #4.b Get present cov set and static variables (don't change) for predictions
      cov <- stack(list.files(file.path(temp, "pres"), full.names = T))
      sta <- cov[[-which(grepl(paste0(c("landuse", "bio"), collapse = "|"), names(cov)))]]
      area_maxSSS <- numeric()
      
      for (k in 1:length(treatments)) {
        
        for (j in 1:length(scens)) { #scens are a combination of RCP * Quartile
          
          #4.c Load data to predict to for different treatments and scenarios
          
          #i Present prediction (current land use and covariates)
          if (treatments[k] == "pre") {
            dyn <- cov[[which(grepl("bio", names(cov)))]]
            luc <- cov[[which(grepl("landuse", names(cov)))]]
          
            #ii 'Aggregated' (Future predictions of land use and climate predictors)
          } else if(treatments[k] == "agg") {
            dyn <- stack(list.files(file.path(temp, paste0("bio", scens[j], "_", country_abbr_preds)), full.names = T))
            luc <- stack(list.files(file.path(temp, paste0("landuse", scens[j], "_", country_abbr_preds)), full.names = T))
            
            #iii Direct: Land use stays the same (current), only climate predictors change
          } else if(treatments[k] == "dir") {
            dyn <- stack(list.files(file.path(temp, paste0("bio", scens[j], "_", country_abbr_preds)), full.names = T))
            luc <- cov[[which(grepl("landuse", names(cov)))]]
            
            #iv Indirect: Land use changes, cliamte predictors stay the same
          } else if(treatments[k] == "ind") {
            dyn <- cov[[which(grepl("bio", names(cov)))]]
            luc <- stack(list.files(file.path(temp, paste0("landuse", scens[j], "_", country_abbr_preds)), full.names = T))
          }
          
          #v Combine the data we just loaded
          cov.new <- stack(sta, dyn, luc)
          cov.new <- stack(cov.new[[which(names(cov.new)%in%names.new)]])
          
          #4.d Predict map
          map <- dismo::predict(mod.final, cov.new)
          rm(cov.new)
          
          #i Apply bioregion constraint for aus (see Methods)
          if (country_abbr == "aus"){
            map[which(!bio_region[]%in%na.omit(unique(bio_region[inds])))] <- NA
          }
          
          #ii Get number of cells above or equal to maxSSS threshold: This is deemed suitable habitat
          area_maxSSS <- c(area_maxSSS, length(which(map[] >= thresh_maxSSS)))
          if(treatments[k] == "pre") { break }
        }
      }
    }
  }
  scens
  rm(dyn, luc, cov.new, mod.final)
  gc()
  list(modresults, area_maxSSS, i)
}
save(results, file = file.path(output_path, paste0("results_", country_abbr, ".RData")))
unlink(temp, recursive = T)
stopCluster(cl)
