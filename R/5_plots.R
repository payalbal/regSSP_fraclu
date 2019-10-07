#required packages
require('matrixStats')
require("pBrackets")
require('RColorBrewer')
require("rasterVis")
require("cividis")
require("sf")

#----------------------------------------#
####---1. LOAD DATA AND GLOAL SETTINGS####
#----------------------------------------#

#1.a File paths and other global parameters
inch <- 0.393701 # one inch figure widths 183 89
path_plots <- file.path(".", "output", "plots")
in_path <- file.path(".", "RData")
out_path <- file.path(".", "output")
fig_path <- file.path(".", "figures")
scenarios <- scens <- c(paste0("q1", c("_26", "_85")), paste0("q2", c("_26", "_85")), paste0("q3", c("_26", "_85")))
treatments <- c("pres", paste0("ind_", scenarios), paste0("dir_", scenarios), paste0("inddir_",scenarios))
ca <- c("aus", "til")

#1.b Get AUC values
auc_list <- list()
exclude_list <- list()
j <- 2
for(j in 1:2){
  res <- readRDS(file.path(out_path, paste0("results_", ca[j], ".rds")))
  auc_list[[j]] <- numeric()
  for(i in 1:length(res)){
    auc_list[[j]][i] <- res[[i]][[1]][,1][which(names(res[[i]][[1]][,1]) == "Test.AUC")]
  }
  exclude_list[[j]] <- which(auc_list[[j]] < 0.7)
}

#1.c Get habitat change results, minus the ones where AUC were insufficient
final_data <- list()
i <- 2
for(i in 1:2){
  res <- readRDS(file.path(out_path, paste0("results_", ca[i], ".rds")))
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
tf <- 0.8
countries <- c("Australia", "Vietnam")

#ii Open device
pdf(file.path(fig_path, "habitat_change.pdf"), height = 10 * inch, width = 18.3 * inch, pointsize = 12)

#iii final plot layout
layout(matrix(c(1, 2, 9, 3, 5,
                1, 2, 9, 3, 6,
                1, 2, 9, 3, 7,
                9, 4, 4, 4, 8), 
              nrow = 4, byrow = T), 
       widths = c(0.1, 0.3,0.05, 0.3,0.5),
       heights = c(0.4, 0.45, 0.45, 0.15))
par(mar = c(0.2,4,0.2,0.2))
par(mgp = c(3, 0.5, 0))

#iv empty plot with axes and some labels
plot(1, ylim = lims, xaxt = "n", type = "n", bty = "n", ylab = "", yaxt = "n")
axis(2, at = c(log(0.125), log(0.25), log(0.5), log(1), log(2), log(4), log(8)), labels = c(expression(paste(over(1,8), " x")),
                                                                                            expression(paste(over(1,4), " x")),
                                                                                            expression(paste(over(1,2), " x")),
                                                                                            expression(paste(1, " x")),
                                                                                            expression(paste(2, " x")),
                                                                                            expression(paste(4, " x")),
                                                                                            expression(paste(8, " x"))),
     las = 2, cex.axis = 0.8, tck = -0.15)

title(ylab = "decrease", cex.lab = 1,
      line = 2.5, adj = 0.3)

title(ylab = "increase", cex.lab = 1,
      line = 2.5, adj = 0.7)


#2.b Plot data points with means, by country (j) and treatment * RCP (i)
for(j in 1:2){
  dp <- final_data[[j]][,which(grepl("q2", colnames(final_data[[j]])))]
  dp <- as.matrix(dp[,c(5, 1, 3, 6, 2, 4)])
  alpha <- 0.1
  means <- colMeans(dp)
  sds <- colSds(dp)
  par(mar = c(0,0,0.2,0))
  plot(1:7, ylim = lims, axes = F, type = "n", xlim = c(0.5, 7.5))
  
  
  for(i in 1:6){
    pts <- dp[,i]
    
    #Only plot points within certain bounds
    pts <- pts[which(pts > log(0.125) & pts < log(8))]
    
    points(rep(x[i], length(pts)) + runif( length(pts), -jit, jit), pts, pch = 20, col = grey(0.7), cex = 1)
    
    #i Means
    segments(x[i] - wth, means[i], x[i] + wth, means[i],
             col= grey(0), border=par("fg"), lty= 1, lwd= 2, xpd=FALSE)
  }
  
  
  title(xlab = countries[j], cex.lab = 1, line = -26)
  of <- 0.1
  brackets(x1 = 0.5, y1 = lims[1] - of, x2 = 3.5, y2 = lims[1]- of, type = 4, ticks = c(-0, -1), h = 0.08)
  brackets(x1 = 4.5, y1 = lims[1] - of, x2 = 7.5, y2 = lims[1]- of, type = 4, ticks = c(0, 1), h = 0.08)
  text(x = 2, y = lims[1] - of - 0.1, "RCP 2.6", cex = tf, xpd = NA)
  text(x = 6, y = lims[1] - of - 0.1, "RCP 8.5", cex = tf, xpd = NA)
  
  
  trmt <- c(1, 2, 3)
  trmt <- c(trmt, trmt)
  posis <- c(1,2,3,5,6,7)
  for (i in 1:6){
    text(x = posis[i], y = log(0.125) - 0.1 - of , trmt[i], cex = tf)
  }
}

par(mar = c(0.5,0,0,0))

plot.new()

#ii Legend
legend("bottom", bty= "n",
       legend = c("1: indirect + direct     2: indirect     3: direct"),
       horiz = T,
       cex = 1,
       border=c(0)
)

#2.c Plot boxplots for AUC
par(mar = c(1.5,8,2,6))
f <- boxplot(auc_list[[1]][-exclude_list[[1]]], auc_list[[2]][-exclude_list[[2]]],
        boxlty = 0,
        staplecol = "black",
        col = "grey",
        boxwex = 0.8,
        boxlty = 0, 
        bty= "n",
        axes = F
)

axis(2, tck = -0.025, at = c(0.7, 0.8, 0.9, 1), labels = c(0.7, 0.8, 0.9, 1))
title(xlab = paste0("n = ", length(auc_list[[1]]) - length(exclude_list[[1]])), adj = 0.1, line = 0)
title(xlab = "Australia", adj = 0.05, line =- 6)
title(xlab = paste0("n = ", length(auc_list[[2]]) - length(exclude_list[[2]])), adj = 0.9, line = 0)
title(xlab = "Vietnam", adj = 0.95, line = - 6)
title(ylab = "AUC", line = 2)


#2.d Plot Varibale contributions
#i get variable contributions
contr_list <- list()
contr_num <- list()
for(j in 1:2){
  preds <- c(readRDS(file.path(out_path, paste0("preds_", ca[j], ".rds"))), "landuse")
  res <- readRDS(file.path(out_path, paste0("results_", ca[j], ".rds")))
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
j <- 1
for(j in 1:2){
  means <- colMeans(contr_list[[j]], na.rm = T)
  ns <- apply(contr_list[[j]], 2, FUN = function (x) {round(length(which(!is.na(x)))/length(x) * 100)})
  means <- means/sum(means) * 100
  names <- names(means)
  names[which(names == "dila")] <- "dist lakes"
  names[which(names == "diri")] <- "dist rivers"
  names(means) <- names
  par(mar = c(0,7,0,3))
  #breaks <- barplot(sort(means), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 25))
  breaks <- barplot(sort(ns), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 100))
  text(sort(ns), breaks, paste0(sort(ns), "%"), adj = 1)
  title(xlab = countries[j], line = -2, adj = 1)
}


par(mar = c(0.1, 7, 0, 3))
par(mgp = c(3, 0.2, 0))
axis(1, xlim = c(0, 100), line = 0.2, tck = -0.025)
plot.new()
title(xlab = "Fraction of models [%]", line = -1.5)
dev.off()

#---------------------#
#####---3. Figure 3####
#---------------------#
countries <- c("Australia", "Vietnam")
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

pdf(file.path(fig_path, "gtap.pdf"), height = 8 * inch, width = 8.9 * inch, pointsize = 12)
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
  title(xlab = countries[i], line = - 1.1)
}
par(mar = c(0,0,0,0))
plot.new()
legend("center", legend = c("Sector output", "Land requirement"), col = c("black", "grey"), lty = c(1, 1), horiz = T, bty = "n", lwd = c(2,2))

dev.off()


#3.b landuse change maps (Also for Supplementary FIgure 2)
#i Global paramerters
cols <- colorRampPalette(c("blue", "yellow", "red"))(101)
breaks <- seq(0, 1, length.out = length(cols))
country <- c("aus", "vnm")
qs <- c("q1", "q2", "q3")
rcps <- c("26", "85")

lu_list <- list()
for(j in 1:2){ #country
  f0 <- readRDS(file.path(in_path, paste0("covariates_", country[j], ".rds")))$landuse
  r_stack <- stack()
  for(k in 1:3){ #quartile
    for(i in 1:2){ #rcp
      name <- paste0("landuse", qs[k], "_", rcps[i], "_", country[j])
      print(name)
      f1 <- readRDS(file.path(out_path, paste0(name, ".rds")))
      dif <- f0 - f1
      dif[which(dif[] != 0 & !is.na(dif[]))] <- 1
      r <- aggregate(dif, fact = 5)
      names(r) <- name
      r_stack <- stack(r_stack, r)
      par(mar = c(0,0,0,0), oma = c(0,0,0,0))
      png(file.path(fig_path, paste0(name, ".png")), bg = "transparent")
      
      print(levelplot(r, col.regions = cols, colorkey = NULL, margin = F, par.settings = list(
        axis.line = list(col = "transparent"),
        strip.background = list(col = 'transparent'), 
        strip.border = list(col = 'transparent')),
        scales = list(draw = FALSE),
        xlab = NULL,
        ylab = NULL,
        at= breaks
      ))
      dev.off()
    }
  }
  lu_list[[j]] <- r_stack
}

#Plot legend image
legend_image <- as.raster(matrix(rev(cols), ncol=1))
pdf(file.path(fig_path, paste0("legend", ".pdf")), height = 4 * inch, width = 2 * inch, pointsize = 12)
par(mar = c(0,0,0,1), oma = c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.1, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex = 0.8)
rasterImage(legend_image, 0, 0, 0.5,1)

dev.off()
#3.c) Land use change per class table
#i Global parameters
rcps <- c("_26", "_85")
countries <- c("vnm", "aus")
mm <- matrix(NA, ncol = 6, nrow = 6)
lu_classes <- c("Urban", "Cropland", "Herbaceous", "Shrubs", "Open forest", "Closed forest")
#ii Open device
pdf(file.path(fig_path, "landuse_table.pdf"), height = 6 * inch, width = 13 * inch, pointsize = 12)

par(mar = c(2, 6, 4, 2))
#iii Get data

for(j in 1:2){
  m <- matrix(NA, ncol = 3, nrow = 6)
  m[, 1] <- as.numeric(readRDS(file.path(".", "output", paste0("final_dmdq2",rcps[j], "_", countries[j], ".rds")))[1,])
  for(i in 1:2){
    m[, i + 1] <- as.numeric(readRDS(file.path(".", "output", paste0("final_dmdq2",rcps[i], "_",  countries[j], ".rds")))[6,])
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
#4.a Get species declining under didfernet scens/treatments (for results section)
countries <- c("aus", "vnm")
i <- 2
mask <- readRDS(file.path(in_path, paste0("mask_", countries[i], ".rds")))

dat <- final_data[[i]]
#iii get species with declines more than 90% under dir and more than 5% under indir effects of RCP8.5 (median predictions)
#i more than 5% decline under indirect imapcts (median predictions)
apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.95))})

#ii more than 90% decline under direct impacts  (median predictions)
apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.05))})

dir_specs <- as.character(dat[,1])[which(exp(dat$dir_q2_85) < 0.05)]
ind_specs <- as.character(dat[,1])[which(exp(dat$ind_q2_85) < 0.95)]

dir_inds <- which(dat$species%in%dir_specs)
ind_inds <- which(dat$species%in%ind_specs)

res <- readRDS(file.path(out_path, paste0("results_", ca[i], ".rds")))
species <- sapply(res, FUN = function (x) {cbind(x[[4]])})
id <- which(species%in%ind_specs)[6]
res[[id]][[4]]
t <- mask
t <- mask - 1
plot(t)
t[res[[id]][[5]]] <- 1
plot(t)
print(species[id])

dir_specs <- as.character(dat[,1])[which(exp(dat$inddir_q2_85) <= 0.5)]
172/nrow(dat)

#4.b Make map of habitat declines across all species
#i Prepare data for plotting
final_maps <- list()
for(i in 1:2){
  country <- ca[i]
  if(country == "til") country <- "vnm"
  mask <- readRDS(file.path(in_path, paste0("mask_", country, ".rds")))
  
  final_maps[[i]] <- stack(mask, mask)
  names(final_maps[[i]]) <- c("indirect", "direct")
  res <- readRDS(file.path(out_path, paste0("results_", ca[i], ".rds")))
  for(j in 1:length(res)){
    if(!res[[j]][[4]]%in%final_data[[i]][,1]) next
    for(k in 5:6){
      if(length(res[[j]][[k]]) != 0) final_maps[[i]][[k-4]][res[[j]][[k]]] <- final_maps[[i]][[k-4]][res[[j]][[k]]] + 1
    }
    print(j)
  }
  final_maps[[i]] <- final_maps[[i]]/nrow(final_data[[i]])
  final_maps[[i]] <- aggregate(final_maps[[i]], fac = 5)
}

#ii Plot individual maps
cols <- colorRampPalette(c("royalblue2", "yellow2", "red2"))(101)
tm <- c('ind', '', "dir")
countries <- c("Australia", "Vietnam")
dev.off()

for(i in 1:2){
  f <- final_maps[[i]] * 100
  for(j in c(1,2)){
    s <- hist(f[[j]])$breaks
    g <- seq(s[1], tail(s,1), length.out = 20)
    cols <- colorRampPalette(c("royalblue2", "yellow2", "red2"))(50)
    png(file.path(fig_path, paste0("habitat_change_", tm[j], "_", ca[i],  ".png")) , bg = "transparent")
    print(levelplot(f[[j]], col.regions = cols, margin = F, colorkey= NULL, par.settings = list(
      axis.line = list(col = "transparent"),
      strip.background = list(col = 'transparent'), 
      strip.border = list(col = 'transparent')),
      scales = list(draw = FALSE),
      xlab = NULL,
      ylab = NULL,
      at = g
    ))
    dev.off()
  }
}

#ii Plot legend for manual insertion in final figure
legend_image <- as.raster(matrix(rev(cols), ncol=1))
pdf(file.path(fig_path, paste0("legend_habitatmaps", ".pdf")), height = 4 * inch, width = 2 * inch, pointsize = 12)
par(mar = c(0,0,0,1), oma = c(0,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.1, y = seq(0,1,l=5), labels = seq(0,1,l=5), cex = 0.8)
rasterImage(legend_image, 0, 0, 0.5,1)
dev.off()

#Land-use Effect sizes
country_abbr <- "vnm"
lu_classes <- c("Urban", "Cropland", "Herbaceous vegetation", "Shrubs", "Open forest", "Closed forest", "Herbaceous wetlands, moss and lichen", "Bare soil and sparse vegetation")
load(file.path(out_path, paste0("landuse_preds", country_abbr, ".RData")))
rm(covariates,ef, obs, pa)

gc()
pdf(file.path(fig_path, paste0("lu_effects_", country_abbr, ".pdf")), height = 18.3 * inch, width = 18.3 * inch, pointsize = 12)
par(mfrow = c(3, 2), mar = par("mar") - c(1, 0.5 , 2, 1))
i <- 1
for(i in 1:(length(lu_classes)-2)){
  sum_mod <- summary(glm.models@models[[i]])$coefficients[-1,]
  suit_preds <- or_names
  cov_names <- paste0("ef_", sprintf("%02d", 1:14))
  xs <- na.omit(match(rownames(sum_mod), cov_names))
  ys <- as.numeric(sum_mod[,1])
  lower <- ys - 1.96 * sum_mod[,2]
  upper <- ys + 1.96 * sum_mod[,2]
  plot(1:14, type = "n", ylim = c(floor(min(lower)), ceiling(max(upper))), xaxt='n', ann = FALSE)
  abline(h = 0, col =alpha("black", alpha = 0.5), lwd = 2)
  arrows(x0 = xs, lower, x1 = xs, upper, length = 0, angle = 90, 
         code = 3, col = alpha("black", alpha = 0.5), lwd = 2)
  points(xs, sum_mod[,1], pch = 16)
  axis(1, c(0,xs), labels = FALSE)
  text(x = c(0,xs), labels = c("", suit_preds[xs]), srt = 45, pos = 1, xpd = TRUE, par("usr")[3] - (par("usr")[4] - par("usr")[3])/10)
  title(ylab="Effect size", line = 2.1)
  title(main=lu_classes[i], line = 1, cex = 0.8)
}
dev.off()

#--------------------------------#
#####---5. Supplementary Plots####
#--------------------------------#
#5.a Global parameters
qs <- c("q1", "q3")
x <- c(1:3, 5:7)
lims <- c(log(0.03), log(8))

#2.a Some parameters
#i GLobal plotting pars
x <- c(1:3, 5:7)
lims <- c(log(0.1), log(8))
padj <- - 0.05
wth <- 0.35
jit <- 0.3
tf <- 0.8
countries <- c("Australia", "Vietnam")
dev.off()
#5.b Supplementary Figure 1
pdf(file.path(fig_path, "sup_habitat_change.pdf"), height = 18.3 * inch, width = 18.3 * inch, pointsize = 12)

layout(matrix(c(1, 2, 8, 3,
                4, 5, 8, 6,
                8, 7, 7, 7), 
              nrow = 3, byrow = T), 
       widths = c(0.1, 0.4, 0.03, 0.4),
       heights = c(0.5, 0.5, 0.05))

for(k in 1:2){
  par(mar = c(0.2,0,0.2,0))
  par(mgp = c(3, 0.5, 0))
  
  plot(1, ylim = lims, xaxt = "n", type = "n", bty = "n", ylab = "", yaxt = "n")
  axis(2, at = c(log(0.125/2), log(0.125), log(0.25), log(0.5), log(1), log(2), log(4), log(8)), 
       labels = c(expression(paste(over(1,16), " x")),
                  expression(paste(over(1,8), " x")),
                  expression(paste(over(1,4), " x")),
                  expression(paste(over(1,2), " x")),
                  expression(paste(1, " x")),
                  expression(paste(2, " x")),
                  expression(paste(4, " x")),
                  expression(paste(8, " x"))),
       las = 2, cex.axis = 0.8, tck = -0.05, line = -5)
  
  title(ylab = "decrease", cex.lab = 1,
        line = -3, adj = 0.3)
  
  title(ylab = "increase", cex.lab = 1,
        line = -3, adj = 0.9)
  
  for(j in 1:2){
    dp <- final_data[[j]][,which(grepl(qs[k], colnames(final_data[[j]])))]
    dp <- as.matrix(dp[,c(5, 1, 3, 6, 2, 4)]) #bring in rioght order
    alpha <- 0.1
    means <- colMeans(dp)
    sds <- colSds(dp)
    par(mar = c(0,0,0.2,0))
    plot(1:7, ylim = lims, axes = F, type = "n", xlim = c(0.5, 7.5))
    
    for(i in 1:6){
      pts <- dp[,i]
      pts <- pts[which(pts > log(0.125/2) & pts < log(8))]
      points(rep(x[i], length(pts)) + runif( length(pts), -jit, jit), pts, pch = 20, col = grey(0.7), cex = 1)
      rect(x[i] - wth, means[i], x[i] + wth, means[i],
           col= grey(0), border=par("fg"), lty= 1, lwd= 3, xpd=FALSE)
    }
    
    tf <- 0.8
    title(xlab = paste0(countries[j], " (", qs[k], ")"), cex.lab = 1, line = -25)
    
    brackets(x1 = 0.5, y1 = lims[1] + 0.1, x2 = 3.5, y2 = lims[1] + 0.1, type = 4, ticks = c(-0, -1), h = 0.08)
    brackets(x1 = 4.5, y1 = lims[1] + 0.1, x2 = 7.5, y2 = lims[1] + 0.1, type = 4, ticks = c(0, 1), h = 0.08)
    text(x = 2, y = lims[1] - 0.1, "RCP 2.6", cex = tf)
    text(x = 6, y = lims[1] - 0.1, "RCP 8.5", cex = tf)
    
    
    trmt <- c(1, 2, 3)
    trmt <- c(trmt, trmt)
    posis <- c(1,2,3,5,6,7)
    for (i in 1:6){
      text(x = posis[i], y = log(0.125/2) - 0.3, trmt[i], cex = tf)
    }
  }
}

par(mar = c(0,0,0,0))

plot.new()

legend("bottom", bty= "n",
       legend = c("1: indirect + direct     2: indirect     3: direct"),
       horiz = T,
       cex = 1,
       border=c(0)
)
dev.off()

#5.b Supplementary FIgure 3 (GTAP output maps for 2071)
s <- read_sf("~/Downloads/TM_WORLD_BORDERS_SIMPL-0.3/", layer = "TM_WORLD_BORDERS_SIMPL-0.3")
t <- read.csv("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/demand calcs/prod and output/qfe_land_all_4oC_2071.csv")
r <- read.csv("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/demand calcs/prod and output/qo_all_4oC_2071.csv")

croptype <- "wht"
borders <- s
borders$QFE85 <- NA
borders$QO85 <- NA
borders$QFE26 <- NA
borders$QO26 <- NA
names(t)[which(!names(t[1,])%in%tolower(borders$ISO3))]
tolower(borders$NAME)[which(!tolower(borders$ISO3)%in%names(t))]

tinds <- na.omit(match(tolower(borders$ISO3) , names(t)))
rinds <- na.omit(match(tolower(borders$ISO3) , names(r)))
finalt <- which(!is.na(match(tolower(borders$ISO3) , names(t))))
finalr <- which(!is.na(match(tolower(borders$ISO3) , names(r))))
borders$QFE85[finalt] <- as.numeric(t[which(t[,1] == croptype), tinds])
borders$QO85[finalr] <- as.numeric(r[which(r[,1] == croptype), rinds])

t <- read.csv("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/demand calcs/prod and output/qfe_land_all_1oC_2071.csv")
r <- read.csv("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/demand calcs/prod and output/qo_all_1oC_2071.csv")

tinds <- na.omit(match(tolower(borders$ISO3) , names(t)))
rinds <- na.omit(match(tolower(borders$ISO3) , names(r)))
finalt <- which(!is.na(match(tolower(borders$ISO3) , names(t))))
finalr <- which(!is.na(match(tolower(borders$ISO3) , names(r))))
borders$QFE26[finalt] <- as.numeric(t[which(t[,1] == croptype), tinds])
borders$QO26[finalr] <- as.numeric(r[which(r[,1] == croptype), rinds])

pdf(file.path(fig_path, "gtap_maps.pdf"))
plot(borders[c("QO26", "QFE26", "QO85", "QFE85")], key.pos = 1, key.length = 0.5, key.width = 0.2)
dev.off()