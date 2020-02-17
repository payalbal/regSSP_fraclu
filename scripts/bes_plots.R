
x <- c('sp', 'raster', 'RColorBrewer')
lapply(x, require, character.only = TRUE)
regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' 
rdata_path <- file.path(regSSP_birds_data, "RData") 
output_path <- file.path(regSSP_birds_data, "output") 

ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))

region = 'vn'



## Barplots for gtap_endowments ####
dmd1 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp1.rds"))
dmd2 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp2.rds"))
dmd3 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp3.rds"))
dmd <- cbind(dmd1[,52], dmd2[,52], dmd3[,52]) # demand in 2070
dmd <- (dmd - mean(dmd))/(max(dmd) - min(dmd)) #mean normalization
dmd <- dmd*100 # convert to percentages
colnames(dmd) <- c("ssp1","ssp2","ssp3")
rownames(dmd) <- c("paddy", "cereal grains", "vegetables,fruits,nuts", "oil seeds", "sugarcane, sugarbeet", "plant-based fibres", "other crops", "cattle,sheep,goats,horses", "forestry")
cols <- brewer.pal(dim(dmd)[1], "RdBu")
par(mar = c(5,4,10,6), oma = c(0,0,0,0))
barplot(dmd, horiz = T, beside = T, col = cols, xlab = "% change in land harvested\n (2019-2070)", ylab = "scenarios")
par(xpd=TRUE)
legend(7,53, rownames(dmd), fill = cols)

## Plot and legend separate
par(mar = c(5,5,2,2), oma = c(0,0,0,0))
barplot(dmd, horiz = T, beside = T, col = cols, xlab = "% change in land harvested (2019-2070)", ylab = "scenarios", cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5)
plot(NULL)
legend(0.2,1, rownames(dmd), fill = cols, cex = 1.2)

## dmd3 only
par(mar = c(5,15,1,1), oma = c(0,0,0,0))
x <- barplot(dmd[,1], horiz = T, beside = T, col = cols, yaxt="n", xlab = "% change in land harvested (2019-2070)", cex.axis = 1.2, cex.names = 1.2, cex.lab = 1.2)
labs <- rownames(dmd)
text(x=-58, y=x, labels = labs, xpd=TRUE, srt=0, cex=1.2, pos=2)


ssp = "ssp3"
temp <- as.data.table(read_xls(file.path(regSSP_birds_data, "gtap_data", "ssps",
                                         paste0(region, "_", ssp, ".xls")), sheet = "land"))
temp <- temp[X__1 %in% gtap_sectors]
temp[, c(2:6, ncol(temp)) := NULL]
temp[, c(20:29) := NULL]
colnames(temp)[1] <- "gtap_sector"
temp <- temp[which(temp$gtap_sector %in% rownames(dmd1)),]
temp <- temp[,c(19)]
temp<- as.matrix(temp)
rownames(temp) <- rownames(dmd1)

dmd.gtap <- temp
dmd.gtap <- cbind(dmd.gtap,temp)

## Total land use change table for vn and aus ####
mm <- matrix(NA, ncol = 6, nrow = 6)
lu_classes <- c("Urban", "Cropland", "Herbaceous", "Shrubs", "Open forest", "Closed forest")
scen_sub <- scens[grep("q2_", scens)][c(2,3)]
regions <- c('vn', 'aus')

for (j in 1:2){
  m <- matrix(NA, ncol = 3, nrow = 6)
  m[, 1] <- as.numeric(readRDS(paste0(output_path, "/finaldmd_", scen_sub[j], "_", regions[j], ".rds"))[1,])
  for(i in 1:2){
    m[, i + 1] <- as.numeric(readRDS(paste0(output_path, "/finaldmd_", scen_sub[i], "_", regions[j], ".rds"))[6,])
  }
  if(j == 1) mm[,1:3] <- m
  if(j == 2) mm[,4:6] <- m
}

## Calculate change relative to entire landscape
rs <- colSums(mm)[c(4,5, 1, 2)]
mn <- cbind(round((mm[,4:6] - mm[,4]), 2), round((mm[,1:3] - mm[,1]), 2))[,-c(1,4)]
mm <- round(t(t(mn)/rs) * 100, 2)
cols1 <- colorRampPalette(c("darkred", "white"))(50)
cols2 <- colorRampPalette(c("white", "darkgreen"))(50)
breaks <- seq(-3, 3, length.out = 101)
new_breaks <- scales::rescale(breaks^3, to = c(-4, 4))
colnames(mm) <- c(c("SSP 1", "SSP 3"), c("SSP 1", "SSP 3"))
rownames(mm) <- lu_classes

## Plot
x <- 1:ncol(mm)
y <- 1:nrow(mm)
centers <- expand.grid(y,x)
par(mar = c(2, 8, 5, 2))
image(x, y, t(mm),
      col = c(cols1, cols2),
      breaks = new_breaks,
      xaxt = 'n', 
      yaxt = 'n',
      xlab = '', 
      ylab = '',
      ylim = c(max(y) + 0.5, min(y) - 0.5)
)

text(centers[,2], centers[,1], c(mm), col= "black", cex = 1.5)
mtext(paste(attributes(mm)$dimnames[[2]]), at=1:ncol(mm), padj = -1, cex = 1.5)
mtext(c("Australia", "Vietnam") , at= c(1.5, 3.5), padj = -3, cex = 1.5)
mtext(attributes(mm)$dimnames[[1]], at=1:nrow(mm), side = 2, las = 1, adj = 1.05, cex = 1.5)


## Landuse maps ####
landuse <- readRDS(paste0(rdata_path, "/covariates_", region, ".rds"))$landuse
par(mfrow=c(1,1), mar = c(0.1, 0.1, 0.1, 0.1))
plot(landuse, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
lu_classes <- c("urban", "crop", "grass", "shrub", "Oforest", "Cforest", "wetland", "bare")
par(xpd = TRUE)
if(region == 'aus'){
  legend(155, -15, legend = lu_classes, fill=ochRe::ochre_palettes$namatjira_qual, bty = "n", cex = 1.2)
} else {
  legend(110, 23, legend = lu_classes, fill=ochRe::ochre_palettes$namatjira_qual, bty = "n", cex = 1.2)
}


lu11 <- readRDS(paste0(output_path, "/predlu_q1_ssp1_", region, ".rds"))
lu21 <- readRDS(paste0(output_path, "/predlu_q2_ssp1_", region, ".rds"))
lu31 <- readRDS(paste0(output_path, "/predlu_q3_ssp1_", region, ".rds"))

lu12 <- readRDS(paste0(output_path, "/predlu_q1_ssp2_", region, ".rds"))
lu22 <- readRDS(paste0(output_path, "/predlu_q2_ssp2_", region, ".rds"))
lu32 <- readRDS(paste0(output_path, "/predlu_q3_ssp2_", region, ".rds"))

lu13 <- readRDS(paste0(output_path, "/predlu_q1_ssp3_", region, ".rds"))
lu23 <- readRDS(paste0(output_path, "/predlu_q2_ssp3_", region, ".rds"))
lu33 <- readRDS(paste0(output_path, "/predlu_q3_ssp3_", region, ".rds"))

plot(stack(lu11,lu21,lu31,lu12,lu22,lu32,lu13,lu23,lu33), legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)

plot(lu11, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
plot(lu31, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
plot(lu13, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)

plot(NULL)
legend("center", legend = lu_classes, fill = ochRe::ochre_palettes$namatjira_qual, bty = "n", cex = 1.7)
# legend("center", legend = lu_classes, fill = ochRe::ochre_palettes$winmar, bty = "n", cex = 1.7)


## Plot difference in land use...
par(mfrow=c(1,1), mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0, 0, 0, 0))
plot(overlay(landuse, lu11, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(landuse, lu31, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(landuse, lu13, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(lu11, lu23, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(NULL)
legend("topright", legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 3)


## Landuse change maps (Also for Supplementary FIgure 2)
#i Global paramerters
cols <- viridis::cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
breaks <- seq(0, 1, length.out = length(cols)) #do we need the sqrt of break points here to?

regions <- c("vn", 'aus')
ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))


lu_list <- list()
# for(j in 1:2){ #regions
for(j in 1) {
  f0 <- readRDS(file.path(rdata_path, paste0("covariates_", regions[j], ".rds")))$landuse
  r_stack <- stack()
  for(k in 1:length(quartiles)){ #quartile
    for(i in 1:length(ssps)){ #rcp
      name <- paste0("predlu_", quartiles[k], "_", ssps[i], "_", regions[j])
      print(name)
      f1 <- readRDS(file.path(output_path, paste0(name, ".rds")))
      dif <- f0 - f1
      dif[which(dif[] != 0 & !is.na(dif[]))] <- 1
      #r <- aggregate(dif, fact = 5)
      r <- focal(dif, w = matrix(1, 3, 3), fun = function(x) {mean(x, na.rm = TRUE)})
      names(r) <- name
      r_stack <- stack(r_stack, r)
      par(mar = c(0,0,0,0), oma = c(0,0,0,0))
      png(file.path(fig_path, paste0(name, ".png")), bg = "transparent")
      #plot(r, col = cols, legend = FALSE, axes=FALSE, box=FALSE)
      
      levelplot(sqrt(r), col.regions = cols, colorkey = NULL, margin = F, par.settings = list(
        axis.line = list(col = "transparent"),
        strip.background = list(col = 'transparent'),
        strip.border = list(col = 'transparent')),
        scales = list(draw = FALSE),
        xlab = NULL,
        ylab = NULL,
        at= breaks
      )
      dev.off()
    }
  }
  lu_list[[j]] <- r_stack
}
