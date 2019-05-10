
#1. Function to download worldclim tiles
get_wctiles <- function(tiles, var, path, ...){
  
  if(missing(path)){
    path <- getwd()
  }
  
  rs <- raster(nrows = 5, ncols = 12, xmn = -180, xmx = 180,
               ymn = -60, ymx = 90)
  rs[] <- 1:length(rs)
  
  tiles_names <- c(paste0(0, 0:11), paste0(1, 0:11), paste0(2, 0:11), paste0(3, 0:11), paste0(4, 0:11))
  
  points <- xyFromCell(rs, which(tiles_names%in%tiles))
  
  message("Loading required worldclim tiles...")
  biotiles <- list()
  for (i in 1:nrow(points)){
    biotiles[[i]] <- getData(name = "worldclim", var = var, path = path, res = 0.5, lon = points[i,1], lat = points[i,2], ...)
  }
  message(paste0("Data successfully downloaded to ", path))
  biotiles
}

#2. Function to merge the downloaded tiles inqto single raster layers
merge_wctiles <- function(biotiles){
  
  if(!is.list(biotiles)) stop("Please provide list with stacks for each tile")
  
  if (length(biotiles) == 1){
    out <- biotiles[[1]]
    message("Only 1 tile detected, merging not necessary")
  }
  
  if(length(biotiles) > 1){
    bionames <- names(biotiles[[1]])
    out <- list()
    message("Merging tiles...")
    for (i in 1:nlayers(biotiles[[1]])){
      b <- biotiles[[1]][[i]]
      message(paste0(substr(bionames[i], 1, nchar(bionames[i])-3)," | " , i, "/", length(bionames)))
      
      for(j in 2:length(biotiles)){
        b <- raster::merge(b, biotiles[[j]][[i]])
      }
      
      out[[i]] <- b
    }
    out <- stack(out)
    names(out) <- bionames
  }
  out
}

#3. Reduce predictor set, always dropping predictor with highest corrleation with anyother predictor
reduce_predset <- function(cors = matrix(), thresh = 0.7) {
  while (min(abs(cors[abs(cors) >= thresh])) != 1){
    values <- cors[which(abs(cors) > thresh)]
    values[values ==1] <- NA
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}

#4. Write raster layers from objects to disk, so maxent can pull them from there rather than having them all loaded into memory
writeToDisk <- function(covariates, folder){
  dir.create(folder, recursive = T)
  writeRaster(covariates, filename = file.path(folder, paste0(names(covariates), ".tif")), bylayer = T, driver = "GTiff", overwrite = T)
}