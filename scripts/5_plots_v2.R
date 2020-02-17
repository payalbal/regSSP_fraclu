#required packages
library('matrixStats')
library("pBrackets")
library('RColorBrewer')
library("rasterVis")
# library("cividis")
library("sf")

#----------------------------------------#
####---1. LOAD DATA AND GLOAL SETTINGS####
#----------------------------------------#

#1.a File paths and other global parameters
inch <- 0.393701 # one inch figure widths 183 89

regSSP_birds_data <- '/Volumes/discovery_data/regSSP_birds_data' 
rdata_path <- file.path(regSSP_birds_data, "RData")
output_path <- file.path(regSSP_birds_data, "output")

fig_path <- file.path(regSSP_birds_data, "figures")
if(!dir.exists(fig_path)){dir.create(fig_path)}


ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))
treatments <- c("pres", paste0("ind_", scens), paste0("dir_", scens), paste0("agg_", scens))

region <- 'til'

# scenarios <- c(paste0("q1", c("_26", "_85")), paste0("q2", c("_26", "_85")), paste0("q3", c("_26", "_85")))
# treatments <- c("pres", paste0("ind_", scenarios), paste0("dir_", scenarios), paste0("inddir_",scenarios))

regions <- c("til", "aus")

#1.b Get AUC values
auc_list <- list()
exclude_list <- list()

# for(j in 1:2){
for(j in 1) {
  res <- readRDS(file.path(output_path, paste0("results2_", regions[j], ".rds")))
  # res <- readRDS("/Volumes/discovery_data/birds_ccimpacts/output/results_til.rds")
  auc_list[[j]] <- numeric()
  for(i in 1:length(res)){
    auc_list[[j]][i] <- res[[i]][[1]][,1][which(names(res[[i]][[1]][,1]) == "Test.AUC")]
  }
  exclude_list[[j]] <- which(auc_list[[j]] < 0.7)
}

#1.c Get habitat change results, minus the ones where AUC were insufficient
final_data <- list()

# for(i in 1:2){
for(i in 1) {
  res <- readRDS(file.path(output_path, paste0("results2_", regions[i], ".rds")))
  species <- sapply(res, FUN = function (x) {cbind(x[[4]])})
  res <- t(sapply(res, FUN = function (x) {cbind(x[[2]])}))
  res <- data.frame(species, res)
  res <- res[-exclude_list[[i]],]
  if(length(which(res[,2] == 0)) != 0) {res <- res[-which(res[,2] == 0),]}
  res[,-1] <- log(res[,-1]/res[,2])
  names(res)[-1] <- treatments
  res <- res[-which(!is.finite(rowSums(res[,-1]))),]
  final_data[[i]] <- res
}
# data <- final_data
# final_data[[1]] <- final_data[[1]][,c(1,2,24,25,26)]
# treatments <- treatments[c(1,23,24,25)]
# final_data

#---------------------#
#####---2. Figure 2####
#---------------------#

#2.a Some parameters
#i GLobal plotting pars
x <- c(1:3, 5:7)
lims <- c(log(0.1), log(8))
padj <- - 0.05
wth <- 0.35
jit <- 0.3
tf <- 1
# countries <- c("Vietnam", "Australia")
countries <- c("Vietnam")

pdf(file.path(fig_path, "habitat_change.pdf"), height = 10 * inch, width = 18.3 * inch, pointsize = 12)
layout(matrix(c(1, 2, 9, 3, 5,
                1, 2, 9, 3, 6,
                1, 2, 9, 3, 7,
                9, 4, 4, 4, 8), 
              nrow = 4, byrow = T), 
       widths = c(0.1, 0.3,0.05, 0.3,0.5),
       heights = c(0.4, 0.45, 0.45, 0.15))
par(mar = c(0.2,4,0.2,0.2))
par(mgp = c(3, 0.5, 0))

plot(1, ylim = lims, xaxt = "n", type = "n", bty = "n", ylab = "", yaxt = "n")
axis(2, at = c(log(0.125), log(0.25), log(0.5), log(1), log(2), log(4), log(8)), labels = c(expression(paste(over(1,8), " x")),
                                                                                            expression(paste(over(1,4), " x")),
                                                                                            expression(paste(over(1,2), " x")),
                                                                                            expression(paste(1, " x")),
                                                                                            expression(paste(2, " x")),
                                                                                            expression(paste(4, " x")),
                                                                                            expression(paste(8, " x"))),
     las = 2, cex.axis = 0.8)

title(ylab = "decrease", cex.lab = 1,
      line = 2.5, adj = 0.3)

title(ylab = "increase", cex.lab = 1,
      line = 2.5, adj = 0.7)


#2.b Plot data points with means, by country (j) and treatment * RCP (i)
# for(j in 1:2){
for(j in 1) {
  dp <- final_data[[j]][,which(grepl("agg_q2", colnames(final_data[[j]])))]
  dp <- as.matrix(dp)
  alpha <- 0.1
  means <- colMeans(dp)
  # sds <- colSds(dp)
  par(mar = c(0,0,0.2,0))
  plot(1:7, ylim = lims, axes = F, type = "n", xlim = c(0.5, 7.5))
  
  
  for(i in 1:3){
    pts <- dp[,i]
    
    #Only plot points within certain bounds
    pts <- pts[which(pts > log(0.125) & pts < log(8))]
    
    points(rep(x[i], length(pts)) + runif( length(pts), -jit, jit), pts, pch = 20, col = grey(0.7), cex = 1)
    
    #i Means
    segments(x[i] - wth, means[i], x[i] + wth, means[i],
             col= grey(0), border=par("fg"), lty= 1, lwd= 2, xpd=FALSE)
  }
  
  
  # title(xlab = countries[j], cex.lab = 1.5, line = -46)
  text(x = 2, y = lims[1] - of - 0.1, "SSP scenario", cex = tf, xpd = NA)

  
  trmt <- c(1, 2, 3)
  trmt <- c(trmt, trmt)
  posis <- c(1,2,3,5,6,7)
  for (i in 1:3){
    text(x = posis[i], y = log(0.125) - 0.1 - of , trmt[i], cex = tf)
  }
}

dev.off()


## Plot boxplots for AUC
par(mar = c(2,8,2,6))
f <- boxplot(auc_list[[1]][-exclude_list[[1]]], #auc_list[[2]][-exclude_list[[2]]],
             boxlty = 0,
             staplecol = "black",
             col = "grey",
             boxwex = 0.8,
             boxlty = 0, 
             bty= "n",
             axes = F
)

axis(2, tck = -0.025, at = c(0.71, 0.8, 0.9, 1), labels = c(0.7, 0.8, 0.9, 1), cex.axis = 1.8)
title(xlab = paste0("n = ", length(auc_list[[1]]) - length(exclude_list[[1]])), adj = 0.1, line = 1, cex.lab = 2)
# title(xlab = "Vietnam", adj = 0.05, line =- 6)
# title(xlab = paste0("n = ", length(auc_list[[2]]) - length(exclude_list[[2]])), adj = 0.9, line = 0)
# title(xlab = "Australia", adj = 0.95, line = - 6)
title(ylab = "AUC", line = 2.5, cex.lab = 2)


## Plot Varibale contributions
contr_list <- list()
contr_num <- list()
# for(j in 1:2){
for(j in 1) {
  preds <- c(readRDS(file.path(output_path, paste0("preds_", regions[j], ".rds"))), "landuse")
  res <- readRDS(file.path(output_path, paste0("results2_", regions[j], ".rds")))
  contr <- matrix(ncol = length(preds), nrow = length(res))
  colnames(contr) <- preds
  
  for (i in 1:length(res)){
    pred_values <- res[[i]][[1]][which(grepl("contribution", rownames(res[[i]][[1]])))]
    pred_names <- rownames(res[[i]][[1]])[which(grepl("contribution", rownames(res[[i]][[1]])))]
    n <- unlist(strsplit(pred_names, split = ".", fixed = T))
    pred_names <- n[seq(1,length(n), 2)]
    contr[i, na.omit(match(pred_names, colnames(contr)))] <- pred_values
  }
  contr_list[[j]] <- contr
  contr_num[[j]] <- apply(contr, 2, FUN = function (x) {length(which(is.na(x)))})/length(res) * 100
}

#ii  Barplots
# for(j in 1:2){
for(j in 1) {
  means <- colMeans(contr_list[[j]], na.rm = T)
  ns <- apply(contr_list[[j]], 2, FUN = function (x) {round(length(which(!is.na(x)))/length(x) * 100)})
  means <- means/sum(means) * 100
  names <- names(means)
  names[which(names == "dila")] <- "dist lakes"
  names[which(names == "diri")] <- "dist rivers"
  names(means) <- names
  names(ns) <- c("slope", "dist to lakes", "dist to rivers", 
                 "prec wettest month", "prec warm quarter", 
                 "prec cold quarter", "mean diurnal range", 
                 "max temp warm month", "temp annual range", "landuse")
  #breaks <- barplot(sort(means), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 25))
  par(mar = c(0,13,0,3))
  breaks <- barplot(sort(ns), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 100), cex.names = 1.5)
  text(sort(ns)-5, breaks, paste0(sort(ns), "%"), adj = 1, cex = 1.5)
  # title(xlab = countries[j], line = -2, adj = 1)
}


#####---3. Figure 3####

regions <- c("Vietnam", "Australia")
#3.a GTAP trajectoires
all <- readRDS(file = file.path(".", "RData", paste0("results_gtap_qo_qfe", ".rds")))
all <- all[c(1:9, 13),]

ncom <- which(rownames(all)%in%c("pdr", "wht", "gro", "v_f", "osd", "c_b", "pfb", "ocr", "ctl", "frs"))[1:8]

#i Define final Layout
ma <- matrix((length(ncom) + 1):(length(ncom) + 1 + 31), #11:50
             nrow = 8, #10
             ncol = 4,
             byrow = T)

ma <- cbind(1:length(ncom), ma) #1:10
ma <- cbind(ma, (max(ma)+1): (max(ma) + nrow(ma)))
ma <- rbind(c((max(ma) + 1):(max(ma) + ncol(ma))) , ma)
ma <- rbind(ma, c((max(ma) + 1):(max(ma) + ncol(ma))))
ma <- rbind(ma, c(rep(max(ma) + 1,3), rep(max(ma) + 2,3)))

ma <- rbind(ma, max(ma) + 1)
ma <- cbind(ma[,c(1:3)], max(ma) + 1, ma[,c(4:6)])
ma[nrow(ma), c(1, 7)] <- max(ma)
ma[nrow(ma),  4] <- max(ma) - 1
ma <- cbind(ma, max(ma))
mw <- 0.3

cols <- c(scales::alpha("black", alpha = 0.05), 
          "grey",
          "black", 
          scales::alpha("black", alpha = 0.5))


#ii Open device, set layout, plot

pdf(file.path(fig_path, "gtap.pdf"), height = 11.5 * inch, width = 8.9 * inch, pointsize = 12)
par(family = "Helvetica")
comnames <- c("paddy rice", "wheat", "cereal grains", "vegetables, fruit, nuts", "oil seeds", "sugar cane, sugar beet", "plant-based fibres", "crops nec", "cattle, sheep, goats, horses", "forestry")[ncom]

layout(ma,
       widths = c(mw,1,1,mw-0.2,1,1,mw + 0.1, mw - 0.1),
       heights = c(mw + 0.1, rep(1, length(ncom)), mw + 0.2, mw + 0.2))

par(mar = c(0.1,0.05,0.1,0.05))

for(i in 1:length(ncom)){
  plot(y = c(-20, 10), x = c(1,1), type = "n", axes = F, xaxt = "n", yaxt = "n", cex = 0.8)
  text(y = -10, x = 1, comnames[i], cex = 0.9, adj = 0, xpd = NA, font = 2)
}

for(i in 1:length(ncom)){
  for(k in 0:1){
    plot(y = c(-20,10), x = c(0, 1), type = "n",axes = F, xaxt = "n", yaxt = "n")
    #rect(0, -20, 1, 10, col = cols[1], lty = 0)
    lines(x = c(0, 1), c(0, all[i, 1+k]), lty = 1, lwd = 2, col = cols[2])
    lines(x = c(0, 1), c(0, all[i, 3+k]), lty = 1, lwd = 2, col = cols[3])
    lines(x = c(0, 1), c(0, 0), lty = 2, lwd = 1, col = cols[4])
  }
  
  for(k in 0:1){
    plot(y = c(-20,10), x = c(0, 1),  type = "n",axes = F)
    #rect(0, -20, 1, 10, col = cols[1], lty = 0)
    lines(x = c(0, 1), c(0, all[i, 5+k]), lty = 1, lwd = 2, col = cols[2])
    lines(x = c(0, 1), c(0, all[i, 7+k]), lty = 1, lwd = 2, col = cols[3])
    lines(x = c(0, 1), c(0, 0), lty = 2, lwd = 1, col = cols[4])
  }
}

par(mar = c(0.1,0.05,0.1,2), mgp = c(2,0.2, 1))
for(i in 1:length(ncom)){
  plot(y = c(-20, 10), x = c(1,1), type = "n", axes = F)
  axis(4, ylim = c(-20, 10), at = c(-20,-10, 0, 10), tcl = -0.25, line = 0.2, col = cols[4], cex.axis = 0.7)
}

par(mar = c(0,0,0,0))
for(i in 1:6){
  plot.new()
  if(i%in%c(2,4)){
    title(xlab = "RCP 2.6", line = -0.9, cex.lab = 1)
  }
  if(i%in%c(3,5)){
    title(xlab = "RCP 8.5", line = -0.9, cex.lab = 1)
  }
}

for(i in 1:6){
  plot.new()
  if(i%in%c(2:5)){
    axis(1, at = c(0,0.2, 0.4, 0.6, 0.8, 1), labels = F, cex.axis = 0.8, tick = T, line = -1, tcl = -0.25, col = cols[4])
    title(xlab = "2018", line = -0.9, adj = 0.05, cex.lab = 1)
    title(xlab = "2070", line = -0.9, adj = 0.95, cex.lab = 1)
  }
}

par(mar = c(0,0,0,0))
for(i in 1:2){
  plot.new()
  title(xlab = regions[i], line = - 1.1)
}
par(mar = c(0,0,0,0))
plot.new()
legend("center", legend = c("Sector output", "Land requirement"), col = c("black", "grey"), lty = c(1, 1), horiz = T, bty = "n", lwd = c(2,2))

dev.off()


## Landuse change maps (Also for Supplementary FIgure 2)
#i Global paramerters
cols <- viridis::cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
breaks <- seq(0, 1, length.out = length(cols)) #do we need the sqrt of break points here to?

regions <- c("vn", 'aus')
ssps <- paste0("ssp", c(1,3))
rcps <- c("45", "85")
quartiles <- c("q2")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))


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
      # png(file.path(fig_path, paste0(name, ".png")), bg = "transparent")
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
      # dev.off()
    }
  }
  lu_list[[j]] <- r_stack
}


#Plot legend image
cols <- viridis::cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
legend_image <- as.raster(matrix(cols, nrow=1))
pdf(file.path(fig_path, paste0("legend", ".pdf")), height = 0.8 * inch, width = 5 * inch, pointsize = 12)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot(c(0,1),c(0,1), type = "n", axes = FALSE,xlab = '', ylab = '')
text(x=sqrt(seq(0,1,l=5)), y = 0.25, labels = seq(0,1,l=5), cex = 0.8)
rasterImage(legend_image, 0, 0.5, 1,1)
dev.off()

## Land use change per class table
rcps <- c("_26", "_85")
regions <- c("vn", "aus")
mm <- matrix(NA, ncol = 6, nrow = 6)
lu_classes <- c("Urban", "Cropland", "Herbaceous", "Shrubs", "Open forest", "Closed forest")

#ii Open device
pdf(file.path(fig_path, "landuse_table.pdf"), height = 6 * inch, width = 13 * inch, pointsize = 12)

par(mar = c(2, 6, 4, 2))
#iii Get data

for(j in 1:2){
  m <- matrix(NA, ncol = 3, nrow = 6)
  m[, 1] <- as.numeric(readRDS(file.path(".", "output", paste0("final_dmdq2",rcps[j], "_", regions[j], ".rds")))[1,])
  for(i in 1:2){
    m[, i + 1] <- as.numeric(readRDS(file.path(".", "output", paste0("final_dmdq2",rcps[i], "_",  regions[j], ".rds")))[6,])
  }
  if(j ==1) mm[,1:3] <- m
  if(j == 2) mm[,4:6] <- m
}

#iv Calculate change relative to entire landscape
rs <- colSums(mm)[c(4,5, 1, 2)]
mn <- cbind(round((mm[,4:6] - mm[,4]), 2), round((mm[,1:3] - mm[,1]), 2))[,-c(1,4)]
mm <- round(t(t(mn)/rs) * 100, 2)
cols1 <- colorRampPalette(c("darkred", "white"))(50)
cols2 <- colorRampPalette(c("white", "darkgreen"))(50)
breaks <- seq(-3, 3, length.out = 101)
new_breaks <- rescale(breaks^3, to = c(-4, 4))
colnames(mm) <- c(c("RCP 2.6", "RCP 8.5"), c("RCP 2.6", "RCP 8.5"))
rownames(mm) <- lu_classes

x <- 1:ncol(mm)
y <- 1:nrow(mm)
centers <- expand.grid(y,x)

image(x, y, t(mm),
      col = c(cols1, cols2),
      breaks = new_breaks,
      xaxt = 'n', 
      yaxt = 'n',
      xlab = '', 
      ylab = '',
      ylim = c(max(y) + 0.5, min(y) - 0.5)
)

text(centers[,2], centers[,1], c(mm), col= "black")

mtext(paste(attributes(mm)$dimnames[[2]]), at=1:ncol(mm), padj = -1)
mtext(c("Australia", "Vietnam") , at= c(1.5, 3.5), padj = -3)
mtext(attributes(mm)$dimnames[[1]], at=1:nrow(mm), side = 2, las = 1, adj = 1.05)

dev.off()


#---------------------#
#####---4. Figure 4####
#---------------------#

## Get species declining under didfernet scens/treatments (for results section)
regions <- c("til", "aus")

reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", regions[i], ".rds")))
dat <- final_data[[i]]
dat <- dat[,which(grepl("agg_q2", colnames(dat)))]
dat <- cbind(final_data[[1]][,1], dat)

## species with declines
## more than 10% decline (median predictions)
apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.9))})
## more than 95% decline (median predictions)
apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.05))})

ssp3_10dec_sp <- as.character(dat[,1])[which(exp(dat$agg_q2_ssp3) < 0.05)]
ssp1_90dec_sp <- as.character(dat[,1])[which(exp(dat$agg_q2_ssp1) < 0.9)]

ssp3_10dec_ind <- which(dat[,1]%in%ssp3_10dec_sp)
ssp1_90dec_ind <- which(dat[,1]%in%ssp1_90dec_sp)


## Make map of habitat declines across all species
res <- readRDS(file.path(output_path, paste0("results2_", regions[i], ".rds")))
spar <- t(sapply(res, "[[", 2))
colnames(spar) <- treatments
rownames(spar) <- sapply(res, "[[", 4)
spar <- spar[,which(grepl("agg_q2_", colnames(spar)))]
spar <- spar[which(rownames(spar) %in% final_data[[1]]$species),]
spar <- na.omit(spar)
spar <- as.data.frame(spar)

# ssps <- paste0("ssp", 1:3)
# rcps <- c("45", "60", "85")
# quartiles <- c("q2", "q1", "q3")
# scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
# treatments <- c("pres", paste0("ind_", scens), paste0("dir_", scens), paste0("agg_", scens))
temp <-  sapply(res, "[", 5)
names(temp) <- sapply(res, "[[", 4)
temp <- temp[which(names(temp)%in%final_data[[1]]$species)]

region = "vn"
final_maps <- list()

if(region == "til") region <- "vn"
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
final_maps <- stack(reg_mask, reg_mask, reg_mask, reg_mask)
names(final_maps) <- c("pres", "ssp1", "ssp2", "ssp3")

for(j in 1:length(temp)){
  if(!names(temp)%in%final_data[[1]]$species) next
  for(k in 1:4){
    layer <- sapply(temp, "[", k)
    if(length(layer[[j]]) != 0) final_maps[[k]][layer[[j]]] <- final_maps[[k]][layer[[j]]] + 1
  }
  print(j)
}
final_maps_raw <- final_maps
final_maps <- final_maps/nrow(final_data[[1]])
final_maps_scaled <- final_maps
for(k in 1:nlayers(final_maps)){
  print(paste0("Resampling output map ", k))
  final_maps[[k]] <- focal(final_maps[[k]], w = matrix(1, 3, 3), fun = function(x) {mean(x, na.rm = TRUE)})
}
names(final_maps) <- c("pres", "ssp1", "ssp2", "ssp3")
    
## Plotting
sum_max <- max(values(final_maps), na.rm = TRUE)
sum_min <- min(values(final_maps), na.rm = TRUE)

par(mar = c(0.2,0,0,0), oma = c(0,0,0,0))
breaks <- seq(0, 0.5, length.out = 20)
cols <- colorRampPalette(c("lightgoldenrod1", "rosybrown2", "royalblue4"))(length(breaks) - 1)
plot(final_maps, legend = F, axes=F, box=F, col = cols)
plot(NULL)
plot(final_maps[[1]], legend.only = T, axes=F, box=F, cex = 1.5, col= cols)


## Plot difference
plot(overlay(final_maps_raw[[2]], final_maps_raw[[1]], fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)



cols <- cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
tm <- c('ind', "dir")
countries <- c("Australia", "Vietnam")
dev.off()
maximum <- sqrt(max(c(stack(final_maps[[1]])[]), stack(final_maps[[2]])[], na.rm = TRUE))


  f <- final_maps

    levelplot(sqrt(f[[1]]), col.regions = cols, margin = F, colorkey= NULL, par.settings = list(
      axis.line = list(col = "transparent"),
      strip.background = list(col = 'transparent'), 
      strip.border = list(col = 'transparent')),
      scales = list(draw = FALSE),
      xlab = NULL,
      ylab = NULL,
      at = breaks
    )


#iii. Legend
cols <- cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
legend_image <- as.raster(matrix(cols, nrow=1))
# pdf(file.path(fig_path, paste0("legend_habmap", ".pdf")), height = 0.8 * inch, width = 5 * inch, pointsize = 12)
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot(c(0,1),c(0,1), type = "n", axes = FALSE,xlab = '', ylab = '')
text(x=sqrt(seq(0,1,l=5)), y = 0.25, labels = seq(0,0.02,l=0.005), cex = 3)
rasterImage(legend_image, 0, 0.5, 1,1)

# for trade analysis
cols <- topo.colors(10)[6:10]
legend_image <- as.raster(matrix(cols, nrow=1))
par(mar = c(0,2,0,12), oma = c(0,0,0,0))
plot(c(0,1),c(0,1), type = "n", axes = FALSE,xlab = '', ylab = '')

text(x=sqrt(seq(0,1,l=5)), y = 0.25, labels = c("0", "0.005", "0.010", "0.015", "0.020"), cex = 2.5, srt =90)
rasterImage(legend_image, 0, 0.5, 1,1)
