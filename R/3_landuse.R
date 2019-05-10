#required packages
library("lulcc")

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

#1.b Get predictors, land use and protected areas
#i Predictors
covariates <- readRDS(files[grepl(paste0("covariates_", country_abbr), files)])

#ii Landuse
lu <- covariates$landuse
obs <- ObsLulcRasterStack(x=lu,categories=c(0:4),
                          labels=c("lu0", "lu1", "lu2", "lu3", "lu4"), t = 0)

#iii Protected areas
pa <- covariates$pa

#iv Drop the ones we want to exclude
covariates <- covariates[[-which(grepl(paste(excludes , collapse = "|"), names(covariates)))]] 



#---------------------------------------#
#####---II. BUILDING SUITABILITY GLMs####
#---------------------------------------#

#2.a) Get values, subset and do correlation analysis.
cov_df <- getValues(covariates)
rows <- nrow(cov_df)
cov_df <- na.omit(cov_df)
subs <- sample(1:nrow(cov_df), 10000)
cors <- cor(cov_df[subs,])
preds <- rownames(reduce_predset(cors))

#2.b) reduce_predset removes predictors from correlated pairs (spearmens > 0.7), see methods in main manuscript.
covariates <- covariates[[which(names(covariates)%in%preds)]]
or_names <- names(covariates)
bio_inds <- which(grepl(or_names, pattern = "bio"))

#2.c) Rename predictors for model
names(covariates) <-  paste0("ef_", sprintf("%02d", 1:(nlayers(covariates))))

#2.d) Partition into training and testing set for GLM - 10% of data used to fit models
part <- partition(x=obs[[1]], size=0.1, spatial=F)
train.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part[["train"]], t=0)

#2.e) Formulate and fit models
cats <- c("lu0","lu1","lu2","lu3","lu4")
forms <- list()
ef <- ExpVarRasterList(covariates)
for(i in 1:5){
  forms[[i]] <- as.formula(paste(cats[i], "~", paste(names(ef), collapse="+")))
}

glm.models <- glmModels(formula=forms, family=binomial, data=train.data, obs=obs)

#2.f) get rid of insignificant predictors and refit
forms2 <- list()
for(i in 1:5){
  first_fit <- summary(glm.models)[[i]]
  names_ef_sig <- rownames(first_fit$coefficients)[which(first_fit$coefficients[,4] <= 0.05)][-1]
  forms2[[i]] <- as.formula(paste(cats[i], "~", paste(names_ef_sig, collapse="+")))
}

glm.models <- glmModels(formula=forms2, family=binomial, data=train.data, obs=obs)

#----------------------------------------#
#####---III. FURTHER MODEL DEFINITIONS####
#----------------------------------------#

#3.a Transition matrix
w <- matrix(1, 5, 5)
w[1, 3] <- 0
w[4, c(1,2,3,5)] <- 0
w[5, c(1,2,3,4,5)] <- 0
w[,5] <- 0

#3.b Scenarios
#i Quartiles under which to predict (main results are under q2, the medians)
qs <- c("q1", "q2", "q3")

#ii two climate change scenarios
scens <- c(26, 85)

#iii Time steps
yrs <- 2070 - 2018
timesteps <- c(12, 22, 32, 42, 52)

#3.c Initialising

#i Global variables
obs2 <- obs
n <- length(which(!is.na(lu[])))
#ii Storage for final demand trajectories (icnluding grass/shrub and forest estimations)
dmd_list <- list()


for (j in 1: length(scens)){
  
  for(k in 1:length(qs)){
    
    #iii Load demand file for current scenario and quartile
    dmd <- as.matrix(readRDS(file.path(".", "output", paste0("landdemand_", country_abbr, "_", scens[j], ".rds"))))
    
    #iv Load required future layers for predictions
    fut <- readRDS(files[grepl(paste0("bio", paste0(c(qs[k], scens[j], country_abbr), collapse = "_")), files)])
    fut <- fut[[which(grepl(paste0(preds, collapse = "|"), names(fut)))]]
    names_fut <- names(covariates)[bio_inds]
    names(fut) <- names_fut
    
    #v Calculate annual increment of dynamic layers (the ones that change, bioclim)
    incr <- (fut - covariates[[bio_inds]])/yrs
    
    #vi Load static (don't change into future) and dynamic layers into separate objects, so they can be recombined into future sets easily 
    cov_stat <- covariates[[-bio_inds]]
    cov_dyn <- covariates[[bio_inds]]
      
    for(i in 1:length(timesteps)){
      
      #vii make timestep specific data stack to predict suitability model to
      covariates_final <- stack(cov_stat, cov_dyn + (incr * timesteps[i]))
      
      #iix determine on the fly demand for grass/shrub and forest land, based on predicted suotability
      cov_df <- na.omit(as.data.frame(covariates_final))
      subs2 <- sample(which(!is.na(covariates_final[[1]][])), 10000)
      pgrass <- mean(predict(glm.models@models[[c(2)]], data = covariates_final[subs2], type = "response"))
      pforest <- mean(predict(glm.models@models[[c(3)]], data = covariates_final[subs2], type = "response"))
      
      sp <- c(pforest, pgrass)
      f <- sp / sum(sp)
      d <- n - sum(dmd[i+1,], na.rm = TRUE)
      
      dmd[i+1,2] <- round(d * f[2])
      dmd[i+1,3] <- round(d * f[1])
      
      #ix reassign original predictor names so they can be found by predict function
      names(covariates_final)[bio_inds] <- names_fut
      ef <- ExpVarRasterList(covariates_final, pattern = "ef")
      
      #3.b Allocate changes in line with demand, suitability model, and allocation order (see Methods)
      
      #i Make model object
      ordered.model <- OrderedModel(obs=obs2,
                                    ef=ef,
                                    models=glm.models,
                                    time=0:1,
                                    demand=dmd[c(i:(i+1)),],
                                    order=c(3, 2 ,0,1,4),
                                    rules = w,
                                    mask = pa)
      
      
      
      #------------------------------#
      #####---IV. ALLOCATE CHANGES####
      #------------------------------#
      ordered_out <- allocate(ordered.model, stochastic = F)
      
      #4.a new land-use map becomes input for next iteration
      obs2 <- ObsLulcRasterStack(ordered_out@output[[2]])
      print(paste0(qs[k], " | Scenario: ", scens[j], " | time step: ", i))
    }
    
    print(table(obs2[]))
    
    #4.b Store outputs
    out <- ordered_out@output[[2]]
    names(out) <- "landuse"
    saveRDS(out, file = paste0(output_path, "/landuse", qs[k], "_", scens[j], "_", country_abbr,  ".rds"))
    saveRDS(dmd, file = paste0(output_path, "/final_dmd", qs[k], "_", scens[j], "_", country_abbr,  ".rds"))
  }
}
