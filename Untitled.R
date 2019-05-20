require(matrixStats)
require("pBrackets")
require(RColorBrewer)
require("rasterVis")

#Main figure
inch <- 0.393701 # one inch figure widths 183 89

path_plots <- file.path(".", "output", "plots")
in_path <- file.path(".", "RData")
out_path <- file.path(".", "output")
fig_path <- file.path(".", "figures")
scenarios <- scens <- c(paste0("q1", c("_26", "_85")), paste0("q2", c("_26", "_85")), paste0("q3", c("_26", "_85")))
treatments <- c("pres", paste0("ind_", scenarios), paste0("dir_", scenarios), paste0("inddir_",scenarios))

ca <- c("aus", "til")

#Get AUC values
auc_list <- list()
exclude_list <- list()
for(j in 1:2){
  res <- readRDS(file.path(out_path, paste0("results_", ca[j], ".rds")))
  auc_list[[j]] <- numeric()
  for(i in 1:length(res)){
    auc_list[[j]][i] <- res[[i]][[1]][,1][which(names(res[[i]][[1]][,1]) == "Test.AUC")]
  }
  exclude_list[[j]] <- which(auc_list[[j]] < 0.7)
}


#Main Figure (habitat chabnge)
final_data <- list()

for(i in 1:2){
  res <- readRDS(file.path(out_path, paste0("results_", ca[i], ".rds")))
  res <- t(sapply(res, FUN = function (x) {cbind(x[[2]])}))
  res <- res[-exclude_list[[i]],]
  if(length(which(res[,1] == 0)) != 0) {res <- res[-which(res[,1] == 0),]}
  diff <- log(res/res[,1])
  colnames(diff) <- treatments
  diff <- diff[-which(!is.finite(rowSums(diff))),]
  final_data[[i]] <- diff
}

lims <- c(log(0.1), log(8))
range(final_data)
x <- c(2:4, 6:8)
cols <- c("#548235", "#C6DFB6", "#70AD47")
dev.off()
layout(matrix(c(1, 2, 9, 3, 5,
                1, 2, 9, 3, 6,
                1, 2, 9, 3, 7,
                9, 4, 4, 4, 8), 
              nrow = 4, byrow = T), 
       widths = c(0.1, 0.3,0.05, 0.3,0.5),
       heights = c(0.4, 0.45, 0.45, 0.1))
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
     las = 2, cex.axis = 0.8, tck = -0.15)

title(ylab = "decrease", cex.lab = 1,
      line = 2.5, adj = 0.3)

title(ylab = "increase", cex.lab = 1,
      line = 2.5, adj = 0.7)

for(j in 1:2){
  dat <- final_data[[j]][,which(grepl("q2", colnames(final_data[[j]])))]
  dat <- dat[,c(5, 1, 3, 6, 2, 4)] #bring in rioght order
  dat <-[,c(1,2,4,5)]
  alpha <- 0.1
  means <- colMeans(dat)
  sds <- colSds(dat)
  countries <- c("Australia", "Vietnam")
  par(mar = c(0,0,0.2,0))
  plot(1:6, 1:6, ylim = lims, axes = F, type = "n")
  rect(0, log(0.5), 9, log(0.125),
       col= scales::alpha("red", alpha = 0.5), border=par("fg"), lty= 1, lwd= 0, xpd=FALSE)
  rect(0, log(1), 9, log(0.5),
       col= scales::alpha("red", alpha = 0.2), border=par("fg"), lty= 1, lwd= 0, xpd=FALSE)
  #rect(0, log(1), 11, log(8),
  #    col= scales::alpha(cols[2], alpha = 0.5), border=par("fg"), lty= 1, lwd= 0, xpd=FALSE)
  wth <- 0.3
  
  for(i in 1:6){
    rect(x[i] - wth, means[i]- 2 * sds[i], x[i] + wth, means[i]+ 2* sds[i],
         col="grey80", border=par("fg"), lty= 0, lwd=par("lwd"), xpd=FALSE)
    rect(x[i] - wth, means[i]- sds[i], x[i] + wth, means[i]+  sds[i],
         col= c(cols,cols)[i], border=par("fg"), lty= 0, lwd=par("lwd"), xpd=FALSE)
    rect(x[i] - wth, means[i], x[i] + wth, means[i],
         col="grey80", border=par("fg"), lty= 1, lwd= 3, xpd=FALSE)
  }
  
  ofs <- 0.4
  text(2-ofs,log(0.9), paste0(round(length(dat[,1][which(dat[,1] < log(1) & dat[,1] >= log(0.5))])/nrow(dat) * 100),"%"), col = "red")
  text(2-ofs, log(0.45), paste0(round(length(dat[,1][which(dat[,1] < log(0.5))])/nrow(dat) * 100),"%"), col = "darkred")
  text(2-ofs, log(1.1), paste0(round(length(dat[,1][which(dat[,1] >= log(1))])/nrow(dat) * 100),"%"), col = "darkgreen")
  
  text(6-ofs, log(0.9), paste0(round(length(dat[,4][which(dat[,4] < log(1) & dat[,4] >= log(0.5))])/nrow(dat) * 100), "%"), col = "red")
  text(6-ofs, log(0.45), paste0(round(length(dat[,4][which(dat[,4] < log(0.5))])/nrow(dat) * 100),"%"), col = "darkred")
  text(6-ofs, log(1.1), paste0(round(length(dat[,4][which(dat[,4] >= log(1))])/nrow(dat) * 100),"%"),col = "darkgreen")
  
  title(xlab = countries[j], cex.lab = 1, line = -26)
  
  brackets(x1 = 1.2, y1 = lims[1], x2 = 4.2, y2 = lims[1], type = 4, ticks = c(-0, -1), h = 0.1)
  brackets(x1 = 5.2, y1 = lims[1], x2 = 8.2, y2 = lims[1], type = 4, ticks = c(0, 1), h = 0.1)
  text(x = 2.7, y = lims[1] - 0.1, "RCP 2.6", cex = 1)
  text(x = 6.7, y = lims[1] - 0.1, "RCP 8.5", cex = 1)
}

par(mar = c(0,0,0,0))

plot.new()
legend("bottom", bty= "n", 
       legend = c("indirect + direct", "indirect", "direct"),
       fill = cols,
       horiz = T, 
       cex = 1,
       border=c(0,0,0)
)


par(mar = c(2,8,2,7))
boxplot(auc_list[[1]][-exclude_list[[1]]], auc_list[[2]][-exclude_list[[2]]],
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
title(xlab = "Australia", adj = 0.05, line =- 5.5)
title(xlab = paste0("n = ", length(auc_list[[2]]) - length(exclude_list[[2]])), adj = 0.9, line = 0)
title(xlab = "Vietnam", adj = 0.95, line = -5.5)
title(ylab = "AUC", line = 2)


#Varibale contributions
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


country <- c("Australia", "Vietnam")
j <- 1
for(j in 1:2){
  means <- colMeans(contr_list[[j]], na.rm = T)
  ns <- apply(contr_list[[j]], 2, FUN = function (x) {round(length(which(!is.na(x)))/length(x) * 100)})
  means <- means/sum(means) * 100
  #names(means)[which(names(means)%in%c("dila", "diri"))] <- c("dist lake", "dist river")
  par(mar = c(0,7,0,3))
  breaks <- barplot(sort(means), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 25))
  text(sort(means), breaks, paste0(ns, "%"), adj = 1)
  title(xlab = country[j], line = -5, adj = 0.9)
}

par(mar = c(0.1, 7, 0, 3))
par(mgp = c(3, 0.2, 0))
axis(1, xlim = c(0, 25), line = 0.2, tck = -0.025)
plot.new()
title(xlab = "% contribution", line = -1)
dev.off()

mask <- readRDS("/Users/simon/Dropbox/PhD/chapter1/birds_ccimpacts/RData/mask_aus.rds")
mask <- mask - 1
pres <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter1/birds_ccimpacts/output/presentation maps/glossy_pre.tif")

fut <- raster("/Users/simon/OneDrive - The University of Melbourne/PhD/chapter1/birds_ccimpacts/output/presentation maps/glossy_agg_q2_85.tif")

fpres <- ffut <- mask
fpres[which(pres[] > 0)] <- 1
ffut[which(fut[] > 0)] <- 1

pdf(paste0(path_plots, "glossy_pres.pdf"))

plot(fpres, axes=FALSE, box=FALSE, legend = FALSE)
dev.off()

pdf(paste0(path_plots, "glossy_fut.pdf"))
plot(ffut, axes=FALSE, box=FALSE, legend = FALSE)
log(length(which(fpres[] > 0))/length(which(ffut[] > 0)))
dev.off()
plot(f1)
plot(f)
plot(f)

summary(final_data[[1]])
df <- as.data.frame(final_data[[2]])
plot(df$dir_q2_85, df$inddir_q2_85, ylim = c(-2, 2), xlim = c(-2, 2))
