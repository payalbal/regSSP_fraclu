#required packages
rm(list = ls())


library("lulcc")
require(raster)
require(glmnet)


#--------------------------------------#
#####---I. LOAD DATA AND DEFINE PARS####
#--------------------------------------#
source(file.path(".", "R", "1_functions.R"))

#1.a File paths to otuputs from preprocessing and some other variables
data_path <- file.path(".", "RData")
output_path <- file.path(".", "output")
country_abbr <- "aus" #aus and vnm
excludes <- c("pa", "popd", "landuse") #(we don't want these in our predictor sets)
files <- list.files(data_path, full.names = T)
mask <- readRDS(file.path(".", "RData", paste0("mask_", country_abbr, ".rds")))


#1.b Get predictors, land use and protected areas
#i Predictors
covariates <- readRDS(files[grepl(paste0("covariates_", country_abbr), files)])


names(covariates)
covariates$landuse[covariates$landuse == 0] <- NA


#ii Landuse
lu_or <- lu <- covariates$landuse
lu_exclude <- c(7,8)
lu[lu%in%lu_exclude] <- NA
obs <- ObsLulcRasterStack(x=lu,categories=c(1:6),
                          labels=paste0("lu", 1:6), t = 0)

#iii Protected areas
pa <- covariates$pa

#mask out classes that cna't change
pa[which(lu_or[]%in%lu_exclude)] <- 0
plot(pa)

#iv Drop the ones we want to exclude
covariates <- covariates[[-which(grepl(paste(excludes , collapse = "|"), names(covariates)))]] 



#---------------------------------------#
#####---II. BUILDING SUITABILITY GLMs####
#---------------------------------------#

#2.a) Get values, subset and do correlation analysis.
cov_df <- getValues(covariates)
rows <- nrow(cov_df)
cov_df <- na.omit(cov_df)
cov_df <- scale(cov_df)
subs <- sample(1:nrow(cov_df), 15000)
cors <- cor(cov_df[subs,])
preds <- rownames(reduce_predset(cors))


#2.b) reduce_predset removes predictors from correlated pairs (spearmens > 0.7), see methods in main manuscript.
covariates <- covariates[[which(names(covariates)%in%preds)]]
or_names <- names(covariates)
bio_inds <- which(grepl(or_names, pattern = "bio"))

cov_cent <- attr(cov_df,c("scaled:center"))[which(colnames(cov_df)%in%or_names)]
cov_scal <- attr(cov_df,c("scaled:scale"))[which(colnames(cov_df)%in%or_names)]

#2.c) Rename predictors for model
names(cov_cent) <- names(cov_scal) <- names(covariates) <-  paste0("ef_", sprintf("%02d", 1:(nlayers(covariates))))
covariates <- stack(scale(covariates, center = cov_cent, scale = cov_scal))

#2.d) Partition into training and testing set for GLM - 10% of data used to fit models
ef <- ExpVarRasterList(covariates)
part <- sample(which(!is.na(obs[[1]][])), 20000)
train.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part, t=0)
train.data <- as.matrix(train.data)


#2.e) Formulate and fit models
cats <- obs@labels
forms <- list()
for(i in 1:length(cats)){
  test <- cv.glmnet(y = train.data[,i], x = train.data[,-c(1:6)], family = "binomial", alpha = 1)
  not_sig <- names(which(coef(test)[,1] == 0))
  if(length(not_sig)==0){
    forms[[i]] <- as.formula(paste(cats[i], "~", paste(names(ef), collapse="+")))
  }else{
    forms[[i]] <- as.formula(paste(cats[i], "~", paste(names(ef)[-which(names(ef)%in%not_sig)], collapse="+")))
  }
  print(i)
}

train.data <- as.data.frame(train.data)
glm.models <- glmModels(formula=forms, family=binomial, data=train.data, obs=obs)
covariates <- readAll(covariates)
save(preds, 
     pa,
     lu,
     obs,
     lu_exclude,
     covariates, 
     lu_or,
     or_names, 
     bio_inds, 
     cov_cent, 
     cov_scal, 
     ef,
     glm.models,
     mask,
     file = paste0("landuse_preds", country_abbr, ".RData")
)


#required packages
rm(list = ls(all = TRUE))
gc()
# system("ps")
# system("pkill -f R")
library("lulcc")
library("doParallel")
#-------------------#
#####---II. MDOEL####
#-------------------#
source(file.path(".", "R", "1_functions.R"))

#1.a File paths to otuputs from preprocessing and some other variables
data_path <- file.path(".", "RData")
output_path <- file.path(".", "output")
country_abbr <- "aus" #aus and vnm
excludes <- c("pa", "popd", "landuse") #(we don't want these in our predictor sets)
files <- list.files(data_path, full.names = T)
mask <- readRDS(file.path(".", "RData", paste0("mask_", country_abbr, ".rds")))

load(paste0("landuse_preds", country_abbr, ".RData"))

#----------------------------------------#
#####---III. FURTHER MODEL DEFINITIONS####
#----------------------------------------#
gc()
#3.a Transition matrix
w <- matrix(1, 6, 6)
w[1,2:6] <- 0
elas = c(1, 0.6,0.6,0.6, 0.7,0.9)

clues.parms <- list(jitter.f=0.002,
                    scale.f=0.000001,
                    max.iter=5000,
                    max.diff= 500,
                    ave.diff= 500)

#3.b Scenarios
#i Quartiles under which to predict (main results are under q2, the medians)
qs <- c("q2", "q1", "q3")
#ii two climate change scenarios
scens <- c(26, 85)


#iii Time steps
yrs <- 2070 - 2018
timesteps <- c(12, 22, 32, 42, 52)

#3.c Initialising
#i Global variables
n <- length(which(!is.na(obs@layers[[1]][])))
logfile <- file.path(output_path, paste0("landuse_", country_abbr, "_log.txt"))
writeLines(c(""), logfile)
#ii Storage for final demand trajectories (icnluding grass/shrub and forest estimations)
dmd_list <- list()
j <- 1
k <- 2
i <- 5
cl <- makeCluster(3)
registerDoParallel(cl)

for (j in 1: length(scens)){
  foreach(k = 1:length(qs), .packages = c("lulcc")) %dopar% {
    #for(k in 1:length(qs)){
    cat(paste("Preparing quantile ", country_abbr, scens[j], k, "\n"), file = logfile, append = T)
    
    dmd <- as.matrix(readRDS(file.path(".", "output", paste0("landdemand_", country_abbr, "_", scens[j], ".rds"))))
    cat('t1', file = logfile, append = T)
    obs2 <- obs
    #iii Load demand file for current scenario and quartile
    #dmd <- matrix(rep(dmd[1,], 6), ncol = 6, byrow = TRUE)
    #iv Load required future layers for predictions
    cat('t2', file = logfile, append = T)
    fut <- readRDS(files[grepl(paste0("bio", paste0(c(qs[k], scens[j], country_abbr), collapse = "_")), files)])
    cat('t3', file = logfile, append = T)
    fut <- fut[[which(grepl(paste0(preds, collapse = "|"), names(fut)))]]
    cat('t4', file = logfile, append = T)
    names_fut <- names(covariates)[bio_inds]
    cat('t5', file = logfile, append = T)
    names(fut) <- names_fut
    cat('t6', file = logfile, append = T)
    fut <- crop(fut, mask)
    cat('t7', file = logfile, append = T)
    fut <- stack(scale(readAll(fut), center = cov_cent[bio_inds], scale = cov_scal[bio_inds]))
    cat('t8', file = logfile, append = T)
    #v Calculate annual increment of dynamic layers (the ones that change, bioclim)
    incr <- (fut - covariates[[bio_inds]])/yrs
    cat('t9', file = logfile, append = T)
    #vi Load static (don't change into future) and dynamic layers into separate objects, so they can be recombined into future sets easily 
    cov_stat <- covariates[[-bio_inds]]
    cat('t10\n', file = logfile, append = T)
    cov_dyn <- covariates[[bio_inds]]
    rm(fut)
    
    
    for(i in 1:length(timesteps)){
      cat(paste("Preparing model ", country_abbr, scens[j], k, timesteps[i], "\n"), file = logfile, append = T)
      cat(paste0(Sys.time(), "\n"),  file = logfile, append = T)
      
      #vii make timestep specific data stack to predict suitability model to
      covariates_final <- stack(cov_stat, cov_dyn + (incr * timesteps[i]))
      gc()
      
      #iix determine on the fly demand for grass, shrub and open and closed forest land, based on predicted suotability
      subs2 <- sample(which(!is.na(covariates_final[[1]][])), 10000)
      pgrass <- mean(predict(glm.models@models[[c(3)]], data = covariates_final[subs2], type = "response"))
      pshrub <- mean(predict(glm.models@models[[c(4)]], data = covariates_final[subs2], type = "response"))
      pforesto <- mean(predict(glm.models@models[[c(5)]], data = covariates_final[subs2], type = "response"))
      pforestc <- mean(predict(glm.models@models[[c(6)]], data = covariates_final[subs2], type = "response"))
      
      sp <- c(pgrass, pshrub, pforesto, pforestc)
      f <- sp / sum(sp)
      d <- n - sum(dmd[i+1,], na.rm = TRUE)
      
      dmd[i+1,3] <- round(d * f[1])
      dmd[i+1,4] <- round(d * f[2])
      dmd[i+1,5] <- round(d * f[3])
      dmd[i+1,6] <- round(d * f[4])
      adj <- rowSums(dmd)[1] - rowSums(dmd)[i + 1]
      col_adj <- sample(1:ncol(dmd), 1)
      dmd[i+1, col_adj] <- dmd[i+1, col_adj] + adj
      if(rowSums(dmd)[1] - rowSums(dmd)[i + 1] != 0) stop("sum of demand not equal to n of cells")
      
      #ix reassign original predictor names so they can be found by predict function
      names(covariates_final)[bio_inds] <- names_fut
      ef <- ExpVarRasterList(readAll(covariates_final), pattern = "ef")
      
      #synch na 9otherwise clues won't converge...
      inds <- unique(c(which(is.na(obs2@layers[[1]][])), which(is.na(ef@maps[[1]][]))))
      for (efs in 1:length(ef)){
        ef@maps[[efs]][inds] <- NA
      }
      pa[inds] <- NA
      obs2@layers[[1]][inds] <- NA
      dmd[i,] <- table(obs2[])
      
      cat(paste("Starting CLUE", country_abbr, scens[j], k, timesteps[i], "\n"), file = logfile, append = T)
      
      clues.model <- CluesModel(obs = obs2, 
                                ef = ef, 
                                models = glm.models, 
                                time = 0:1,
                                demand=dmd[c(i:(i+1)),], 
                                rules = w, 
                                params = clues.parms, 
                                elas = elas,
                                mask = pa)
      st <- system.time(ordered_out <- allocate(clues.model))
      cat(paste(st, "\n"), file = logfile, append = T)
      
      #3.b Allocate changes in line with demand, suitability model, and allocation order (see Methods)
      
      # #i Make model object
      # ordered.model <- OrderedModel(obs=obs2,
      #                               ef=ef,
      #                               models=glm.models,
      #                               time=0:1,
      #                               demand=dmd[c(i:(i+1)),],
      #                               order=c(1:6),
      #                               rules = w,
      #                               mask = pa)
      
      
      #------------------------------#
      #####---IV. ALLOCATE CHANGES####
      #------------------------------#
      #ordered_out <- allocate(ordered.model, stochastic = F)
      
      #4.a new land-use map becomes input for next iteration
      obs2 <- ObsLulcRasterStack(ordered_out@output[[2]])
      
      print(paste0(qs[k], " | Scenario: ", scens[j], " | time step: ", i))
    }
    rm(ef)
    #4.b Store outputs
    out <- ordered_out@output[[2]]
    sum(table(out[]))
    out[which(lu_or[]%in%lu_exclude[1])] <- 7
    out[which(lu_or[]%in%lu_exclude[2])] <- 8
    names(out) <- "landuse"
    saveRDS(readAll(out), file = paste0(output_path, "/landuse", qs[k], "_", scens[j], "_", country_abbr,  ".rds"))
    saveRDS(dmd, file = paste0(output_path, "/final_dmd", qs[k], "_", scens[j], "_", country_abbr,  ".rds"))
  }
}

#test that landuse is fine
country_abbr <- "aus"
x <- list.files(output_path, pattern = country_abbr,  full.names = TRUE)
xr <- x[grepl("landuseq", x)]
xd <- x[grepl("final_dmd", x)]
landuse <- readRDS(paste0("./RData/covariates_", country_abbr, ".rds"))
landuse <- landuse$landuse
lu_tabled <- table(landuse[])
rs <- list()
for(i in 1:(length(x)-2)){
  r <- readRDS(xr[i])
  rs[[i]] <- r
  d <- readRDS(xd[i])
  print(all(table(r[])[1:6] == d[6,]))
  print(d[6,])
}
rs_df <- as.data.frame(stack(rs))
for(j in 1:6){
  for(i in 1:6){
    print(identical(rs_df[,j], rs_df[,i]))
  }
}
plot(obs2@layers[[1]])

