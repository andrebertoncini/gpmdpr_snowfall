#Code to process GPM-DPR data (version 2A.GPM.DPR/V9) and generate spatial estimates

setwd("Path to raw netcdf GPM-DPR files")

library(ncdf4)
library(raster)
library(RColorBrewer)
library(hdf5r)


#Id of events with precipitation in the study area

id <- c(39, 50, 55, 69, 78, 90, 94, 97, 116, 118, 126, 129, 132, 133, 134, 158, 164, 172, 189, 200, 209, 213, 219, 236, 253, 254, 257, 260, 261, 270, 298, 300, 307, 311)


#Load input files

bin_free <- read.csv("CSV file with GPM clutter free bin information")[id,"add col number"]

mrr_W_input <- read.csv("CSV file with MRR Doppler velocity in the 1000-m layer above clutter free bin")

mrr_vel <- mrr_W_input[id, 2]*100

provinces <- shapefile("Shapefile with province boundaries")


files <- list.files(getwd())

ka_stack <- stack()


#Zm-Ka maps

#Load GPM-DPR netcdf files

for (k in 1:34) {
  
  netcdf <- nc_open(files[id[k]])
  
  attributes(netcdf$var)
  
  
  long <- as.vector(ncvar_get(netcdf, attributes(netcdf$var)$names[29]))
  lat <- as.vector(ncvar_get(netcdf, attributes(netcdf$var)$names[30]))
  
  
  #Generate spatial clutter free bin
  
  binClutterFreeBottom <- as.vector(ncvar_get(netcdf, attributes(netcdf$var)$names[3]))
  
  binClutterFreeBottom_df <- data.frame(long, lat, binClutterFreeBottom)
  
  binClutterFreeBottom_aoi <- binClutterFreeBottom_df[which(binClutterFreeBottom_df$long > -117 & binClutterFreeBottom_df$long < -114 & binClutterFreeBottom_df$lat > 50 & binClutterFreeBottom_df$lat < 52),]
  
  colnames(binClutterFreeBottom_aoi) <- c("x", "y", "z")
  
  e <- extent(binClutterFreeBottom_aoi[,1:2])
  
  r <- raster(res = 0.07, xmn = -117, xmx = -114, ymn = 50, ymx = 52)
  
  binClutterFreeBottom_raster <- rasterize(binClutterFreeBottom_aoi[, 1:2], r, binClutterFreeBottom_aoi[,3], fun=mean)
  
  plot(binClutterFreeBottom_raster)
  
  
  rm("binClutterFreeBottom")
  rm("binClutterFreeBottom_df")
  
  
  #Generate spatial precipitation occurrence flag
  
  flagPrecip <- as.vector(ncvar_get(netcdf, attributes(netcdf$var)$names[6]))
  
  flagPrecip_df <- data.frame(long, lat, flagPrecip)
  
  flagPrecip_aoi <- flagPrecip_df[which(flagPrecip_df$long > -117 & flagPrecip_df$long < -114 & flagPrecip_df$lat > 50 & flagPrecip_df$lat < 52),]
  
  colnames(flagPrecip_aoi) <- c("x", "y", "z")
  
  e <- extent(flagPrecip_aoi[,1:2])
  
  r <- raster(res = 0.07, xmn = -117, xmx = -114, ymn = 50, ymx = 52)
  
  flagPrecip_raster <- rasterize(flagPrecip_aoi[, 1:2], r, flagPrecip_aoi[,3], fun=mean)
  
  plot(flagPrecip_raster > 0)
  points(-115.213909374496, 50.8259838210696)
  
  rm("flagPrecip")
  rm("flagPrecip_df")
  
  
  #Generate spatial temperature fields in the 1000-m layer above clutter free bin
  
  airTemperature <- ncvar_get(netcdf, attributes(netcdf$var)$names[1])
  
  airTemperature_matrix <- matrix(nrow = 8, ncol = length(binClutterFreeBottom_aoi$z))
  
  
  for (j in (bin_free[k] - 7):bin_free[k]) {
    
    airTemperature_vec <- as.vector(airTemperature[j,,])
    
    airTemperature_df <- data.frame(long, lat, airTemperature_vec)
    
    airTemperature_aoi <- airTemperature_df[which(airTemperature_df$long > -117 & airTemperature_df$long < -114 & airTemperature_df$lat > 50 & airTemperature_df$lat < 52),]
    
    airTemperature_matrix[j-(bin_free[k] - 8),] <- airTemperature_aoi$airTemperature_vec
    
  }
  
  
  airTemperature_aoi <- data.frame(airTemperature_aoi$long, airTemperature_aoi$lat, colMeans(airTemperature_matrix))
  
  colnames(airTemperature_aoi) <- c("x", "y", "z")
  
  e <- extent(airTemperature_aoi[,1:2])
  
  r <- raster(res = 0.07, xmn = -117, xmx = -114, ymn = 50, ymx = 52)
  
  airTemperature_aoi_raster <- rasterize(airTemperature_aoi[, 1:2], r, airTemperature_aoi[,3], fun=mean) - 273.15
  
  plot(airTemperature_aoi_raster)
  
  rm("airTemperature")
  
  
  #Generate spatial measured Zm-Ka fields in the 1000-m layer above clutter free bin
  
  zFactorMeasured <- ncvar_get(netcdf, attributes(netcdf$var)$names[4])
  
  zFactorMeasured_Ka_matrix <- matrix(nrow = 8, ncol = length(binClutterFreeBottom_aoi$z))
  
  for (j in (bin_free[k] - 7):bin_free[k]) {
    
    zFactorMeasured_Ka_vec <- as.vector(zFactorMeasured[2,j,,])
    
    zFactorMeasured_Ka_df <- data.frame(long, lat, zFactorMeasured_Ka_vec)
    
    zFactorMeasured_Ka_aoi <- zFactorMeasured_Ka_df[which(zFactorMeasured_Ka_df$long > -117 & zFactorMeasured_Ka_df$long < -114 & zFactorMeasured_Ka_df$lat > 50 & zFactorMeasured_Ka_df$lat < 52),]
    
    zFactorMeasured_Ka_matrix[j-(bin_free[k] - 8),] <-  zFactorMeasured_Ka_aoi$zFactorMeasured_Ka_vec
    
  }
  
  
  zFactorMeasured_Ka_matrix <- ifelse(zFactorMeasured_Ka_matrix < 0, NA, zFactorMeasured_Ka_matrix)
  
  zFactorMeasured_Ka_aoi <- data.frame(zFactorMeasured_Ka_aoi$long, zFactorMeasured_Ka_aoi$lat, colMeans(zFactorMeasured_Ka_matrix, na.rm = T))
  
  colnames(zFactorMeasured_Ka_aoi) <- c("x", "y", "z")
  
  e <- extent(zFactorMeasured_Ka_aoi[,1:2])
  
  r <- raster(res = 0.07, xmn = -117, xmx = -114, ymn = 50, ymx = 52)
  
  zFactorMeasured_Ka_aoi_raster <- rasterize(zFactorMeasured_Ka_aoi[, 1:2], r, zFactorMeasured_Ka_aoi[,3], fun=mean)
  
  plot(zFactorMeasured_Ka_aoi_raster)
  
  rm("zFactorMeasured")
  
  
  zFactorMeasured_Ka_aoi_raster[flagPrecip_raster == 0] <- NA
  zFactorMeasured_Ka_aoi_raster[zFactorMeasured_Ka_aoi_raster < 12] <- NA
  plot(zFactorMeasured_Ka_aoi_raster)
  
  
  #Generate spatial snowfall fields in the 1000-m layer above clutter free bin using Zm-Ka
  
  airTemperature_aoi_raster[airTemperature_aoi_raster > 0] <- NA
  
  N0_aoi <- (5.65*(10^5)*exp((-0.107*airTemperature_aoi_raster)))/100000000
  
  plot(N0_aoi)
  
  a_coef_aoi <- airTemperature_aoi_raster
  
  a_coef_aoi[airTemperature_aoi_raster < -5] <- ((3.36*(10000/((N0_aoi[airTemperature_aoi_raster < -5]^(2/3))*(mrr_vel[k]^(5/3)))))*0.189)
  
  a_coef_aoi[airTemperature_aoi_raster >= -5] <- ((5.32*(10000/((N0_aoi[airTemperature_aoi_raster >= -5]^(2/3))*(mrr_vel[k]^(5/3)))))*0.189)
  
  b_exp <- 1.67
  
  plot(a_coef_aoi)
  

  precip_Ka_aoi <- ((10^(zFactorMeasured_Ka_aoi_raster/10))/a_coef_aoi)^(1/b_exp)
  
  ka_stack <- stack(ka_stack, precip_Ka_aoi)
  
}


#Generate CORRA spatial estimates

setwd("Path to CORRA HDF5 files")


files <- list.files()

corra_stack <- stack()

for (j in 1:34) {
  
  hdf5_file <- H5File$new(files[id[j]], mode = "r")
  
  dataset = list.datasets(hdf5_file)
  
  
  precipTotRate <- readDataSet(hdf5_file[["/KuKaGMI/precipTotRate"]])
  
  long <- readDataSet(hdf5_file[["/KuKaGMI/Longitude"]])
  
  lat <- readDataSet(hdf5_file[["/KuKaGMI/Latitude"]])
  
  
  precip_bin_free <- as.vector(precipTotRate[floor(bin_free[j]/2),,])[which(lat > 50 & lat < 52 & long < - 114 & long > -117)]
  
  precip_bin_free <- ifelse(precip_bin_free <= 0, NA, precip_bin_free)
  
  
  precip_bin_free_m1 <- as.vector(precipTotRate[(floor(bin_free[j]/2) - 1),,])[which(lat > 50 & lat < 52 & long < - 114 & long > -117)]
  
  precip_bin_free_m1 <- ifelse(precip_bin_free_m1 <= 0, NA, precip_bin_free_m1)
  
  
  precip_bin_free_m2 <- as.vector(precipTotRate[(floor(bin_free[j]/2) - 2),,])[which(lat > 50 & lat < 52 & long < - 114 & long > -117)]
  
  precip_bin_free_m2 <- ifelse(precip_bin_free_m2 <= 0, NA, precip_bin_free_m2)
  
  
  precip_bin_free_m3 <- as.vector(precipTotRate[(floor(bin_free[j]/2) - 3),,])[which(lat > 50 & lat < 52 & long < - 114 & long > -117)]
  
  precip_bin_free_m3 <- ifelse(precip_bin_free_m3 <= 0, NA, precip_bin_free_m3)
  
  
  precip_accum_layer <- (precip_bin_free + precip_bin_free_m1 + precip_bin_free_m2 + precip_bin_free_m3)/4
  
  
  precipTotRate_df <- data.frame(as.vector(long)[which(lat > 50 & lat < 52 & long < - 114 & long > -117)],
                                 as.vector(lat)[which(lat > 50 & lat < 52 & long < - 114 & long > -117)],
                                 precip_accum_layer)
  
  
  colnames(precipTotRate_df) <- c("x", "y", "z")
  
  
  e <- extent(c(-117, -114, 50, 52))
  
  r <- raster(e, res = 0.07)
  
  precipTotRate_raster <- rasterize(precipTotRate_df[, 1:2], r, precipTotRate_df[,3], fun=mean)
  
  
  corra_stack <- stack(corra_stack, precipTotRate_raster)

  
  rm("precipTotRate")
  rm("hdf5_file")
  
}



#IMERG spatial estimates

setwd("Path to IMERG GeoTiff files")


files_imerg <- list.files(getwd())

imerg_stack <- stack()


for (i in 1:34) {
  
  imerg_raster <- raster(files_imerg[i])
  
  imerg_raster[imerg_raster == 0] <- NA
  
  imerg_stack <- stack(imerg_stack, imerg_raster)
  
}



#Accumulated snowfall maps

rockies <- shapefile("Shapefile with Canadian Rockies boundaries")

palette_1 <- brewer.pal(9, "YlGnBu")

ka_stack_masked <- mask(ka_stack, rockies)

corra_stack_masked <- mask(corra_stack, rockies)

imerg_stack_masked <- mask(imerg_stack, rockies)


ka_accum <- calc(ka_stack_masked, sum, na.rm = T)

ka_accum[ka_accum == 0] <- NA


#Zm-Ka accumulated snowfall map

tiff(paste0("Path to save map.tiff"), width = 5, height = 5, units = "in", res = 500)

par(cex = 1.5)

plot(ka_accum, zlim = c(0, 29), col = palette_1)
plot(provinces, add = T, bg = "transparent", border = "black")
plot(rockies, add = T, bg = "transparent", border = "purple")
points(-115.213909374496, 50.8259838210696)
text(-115.213909374496, 50.8259838210696, "Fortress", pos = 4, col = "black", cex = 1.5)

dev.off()


#CORRA accumulated snowfall map

corra_accum <- calc(corra_stack_masked, sum, na.rm = T)

corra_accum[corra_accum == 0] <- NA


tiff(paste0("Path to save map.tiff"), width = 5, height = 5, units = "in", res = 500)

par(cex = 1.5)

plot(corra_accum, zlim = c(0, 29), col = palette_1)
plot(provinces, add = T, bg = "transparent", border = "black")
plot(rockies, add = T, bg = "transparent", border = "purple")
points(-115.213909374496, 50.8259838210696)
text(-115.213909374496, 50.8259838210696, "Fortress", pos = 4, col = "black", cex = 1.5)

dev.off()


#IMERG accumulated snowfall map

imerg_accum <- calc(imerg_stack_masked, sum, na.rm = T)

imerg_accum[imerg_accum == 0] <- NA


tiff(paste0("Path to save map.tiff"), width = 5, height = 5, units = "in", res = 500)

par(cex = 1.5)

plot(imerg_accum, zlim = c(0, 29), col = palette_1)
plot(provinces, add = T, bg = "transparent", border = "black")
plot(rockies, add = T, bg = "transparent", border = "purple")
points(-115.213909374496, 50.8259838210696)
text(-115.213909374496, 50.8259838210696, "Fortress", pos = 4, col = "black", cex = 1.5)

dev.off()


#CDF plots using pixel values of each snowfall product

ka_stack_vec <- as.vector(ka_stack_masked)

corra_stack_vec <- as.vector(corra_stack_masked)

imerg_stack_vec <- as.vector(imerg_stack_masked)


tiff(paste0("Path to save plot.tiff"), width = 5, height = 5, units = "in", res = 500)

plot(ecdf(ka_stack_vec), col = "blue", xlab = "Snowfall Rate [mm/h]", ylab = "Cumulative Distribution Function [ ]", main = " ", lwd = 1.5, xlim = c(0,11))
plot(ecdf(corra_stack_vec), add = T, col = "green", do.points = F, lwd = 1.5)
plot(ecdf(imerg_stack_vec), add = T, col = "purple", do.points = F, lwd = 1.5)

legend(7, 0.3, legend = c("Zm-Ka", "CORRA", "IMERG"), lty = c(1,1,1), col = c("blue", "green", "purple"))

dev.off()

