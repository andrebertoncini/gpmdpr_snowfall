#Code to process GPM-DPR data (version 2A.GPM.DPR/V9) and extract data at the station coordinate
#By Andre Bertoncini
#Centre for Hydrology - University Saskatchewan


setwd("Set working directory with GPM-DPR netcdf files")


library(ncdf4)

files <- list.files(getwd())


#Load GPM-DPR netcdf files and extract variable profiles at station coordinates

for (j in 1:315) {
  
  netcdf <- nc_open(files[j])
  
  attributes(netcdf$var)
  
  
  long <- ncvar_get(netcdf, attributes(netcdf$var)$names[29])
  lat <- ncvar_get(netcdf, attributes(netcdf$var)$names[30])
  
  
  grid_point <- which.min(sqrt(abs(long - round(-115.198195978996,4))^2 + abs(lat - round(50.8259838210696,5))^2))
  
  
  output_matrix <- matrix(ncol = 32, nrow = 176)
  
  
  airTemperature <- ncvar_get(netcdf, attributes(netcdf$var)$names[1])
  airTemperature_prof <- vector()
  
  for (i in 1:176) {
    
    airTemperature_prof[i] <- ((as.vector(airTemperature[i,,])[grid_point]) - 273.15)
    
  }
  
  rm("airTemperature")
  output_matrix[,1] <- airTemperature_prof
  
  
  binZeroDeg <- ncvar_get(netcdf, attributes(netcdf$var)$names[2])
  binZeroDeg_surf <- binZeroDeg[grid_point]
  rm("binZeroDeg")
  output_matrix[176,2] <- binZeroDeg_surf
  
  
  binClutterFreeBottom <- ncvar_get(netcdf, attributes(netcdf$var)$names[3])
  binClutterFreeBottom_surf <- binClutterFreeBottom[grid_point]
  rm("binClutterFreeBottom")
  output_matrix[176,3] <- binClutterFreeBottom_surf
  
  
  zFactorMeasured <- ncvar_get(netcdf, attributes(netcdf$var)$names[4])
  zFactorMeasured_prof_ka <- vector()
  
  for (i in 1:176) {
    
    zFactorMeasured_prof_ka[i] <- as.vector(zFactorMeasured[2,i,,])[grid_point]
    
  }
  
  output_matrix[,4] <- zFactorMeasured_prof_ka
  
  
  zFactorMeasured_prof_ku <- vector()
  
  for (i in 1:176) {
    
    zFactorMeasured_prof_ku[i] <- as.vector(zFactorMeasured[1,i,,])[grid_point]
    
  }
  
  rm("zFactorMeasured")
  output_matrix[,31] <- zFactorMeasured_prof_ku
  
  
  flagPrecip <- ncvar_get(netcdf, attributes(netcdf$var)$names[6])
  flagPrecip_surf <- flagPrecip[grid_point]
  rm("flagPrecip")
  output_matrix[176,5] <- flagPrecip_surf
  
  
  flagBB <- ncvar_get(netcdf, attributes(netcdf$var)$names[7])
  flagBB_surf <- flagBB[grid_point]
  rm("flagBB")
  output_matrix[176,6] <- flagBB_surf
  
  
  binBBBottom <- ncvar_get(netcdf, attributes(netcdf$var)$names[8])
  binBBBottom_surf <- binBBBottom[grid_point]
  rm("binBBBottom")
  output_matrix[176,7] <- binBBBottom_surf
  
  
  binBBTop <- ncvar_get(netcdf, attributes(netcdf$var)$names[9])
  binBBTop_surf <- binBBTop[grid_point]
  rm("binBBTop")
  output_matrix[176,8] <- binBBTop_surf
  
  
  binBBPeak <- ncvar_get(netcdf, attributes(netcdf$var)$names[10])
  binBBPeak_surf <- binBBPeak[grid_point]
  rm("binBBPeak")
  output_matrix[176,9] <- binBBPeak_surf
  
  
  precipRateNearSurface <- ncvar_get(netcdf, attributes(netcdf$var)$names[11])
  precipRateNearSurface_surf <- precipRateNearSurface[grid_point]
  rm("precipRateNearSurface")
  output_matrix[176,10] <- precipRateNearSurface_surf
  
  
  precipRate <- ncvar_get(netcdf, attributes(netcdf$var)$names[12])
  precipRate_prof <- vector()
  
  for (i in 1:176) {
    
    precipRate_prof[i] <- as.vector(precipRate[i,,])[grid_point]
    
  }
  
  rm("precipRate")
  output_matrix[,11] <- precipRate_prof
  
  
  zFactorFinal <- ncvar_get(netcdf, attributes(netcdf$var)$names[13])
  zFactorFinal_prof_ka <- vector()
  
  for (i in 1:176) {
    
    zFactorFinal_prof_ka[i] <- as.vector(zFactorFinal[2,i,,])[grid_point]
    
  }
  
  output_matrix[,12] <- zFactorFinal_prof_ka
  
  
  zFactorFinal_prof_ku <- vector()
  
  for (i in 1:176) {
    
    zFactorFinal_prof_ku[i] <- as.vector(zFactorFinal[1,i,,])[grid_point]
    
  }
  
  rm("zFactorFinal")
  output_matrix[,32] <- zFactorFinal_prof_ku
  
  
  paramDSD <- ncvar_get(netcdf, attributes(netcdf$var)$names[14])
  paramDSD_prof <- vector()
  
  for (i in 1:176) {
    
    paramDSD_prof[i] <- as.vector(paramDSD[2,i,,])[grid_point]
    
  }
  
  rm("paramDSD")
  output_matrix[,13] <- paramDSD_prof
  
  
  zFactorFinalESurface <- ncvar_get(netcdf, attributes(netcdf$var)$names[15])
  zFactorFinalESurface_surf <- zFactorFinalESurface[2,,][grid_point]
  rm("zFactorFinalESurface")
  output_matrix[176,14] <- zFactorFinalESurface_surf
  
  
  precipRateAve24 <- ncvar_get(netcdf, attributes(netcdf$var)$names[16])
  precipRateAve24_surf <- precipRateAve24[grid_point]
  rm("precipRateAve24")
  output_matrix[176,15] <- precipRateAve24_surf
  
  
  zFactorFinalNearSurface <- ncvar_get(netcdf, attributes(netcdf$var)$names[17])
  zFactorFinalNearSurface_surf <- zFactorFinalNearSurface[2,,][grid_point]
  rm("zFactorFinalNearSurface")
  output_matrix[176,16] <- zFactorFinalNearSurface_surf
  
  
  precipRateESurface <- ncvar_get(netcdf, attributes(netcdf$var)$names[18])
  precipRateESurface_surf <- precipRateESurface[grid_point]
  rm("precipRateESurface")
  output_matrix[176,17] <- precipRateESurface_surf
  
  
  phaseNearSurface <- ncvar_get(netcdf, attributes(netcdf$var)$names[19])
  phaseNearSurface_surf <- phaseNearSurface[grid_point]
  rm("phaseNearSurface")
  output_matrix[176,18] <- phaseNearSurface_surf
  
  
  ScanTime_Minute <- ncvar_get(netcdf, attributes(netcdf$var)$names[20])
  ScanTime_Minute_surf <- ScanTime_Minute[round(grid_point/49)]
  rm("ScanTime_Minute")
  output_matrix[176,19] <- ScanTime_Minute_surf
  
  
  ScanTime_Year <- ncvar_get(netcdf, attributes(netcdf$var)$names[21])
  ScanTime_Year_surf <- ScanTime_Year[round(grid_point/49)]
  rm("ScanTime_Year")
  output_matrix[176,20] <- ScanTime_Year_surf
  
  
  ScanTime_Second <- ncvar_get(netcdf, attributes(netcdf$var)$names[22])
  ScanTime_Second_surf <- ScanTime_Second[round(grid_point/49)]
  rm("ScanTime_Second")
  output_matrix[176,21] <- ScanTime_Second_surf
  
  
  ScanTime_Hour <- ncvar_get(netcdf, attributes(netcdf$var)$names[23])
  ScanTime_Hour_surf <- ScanTime_Hour[round(grid_point/49)]
  rm("ScanTime_Hour")
  output_matrix[176,22] <- ScanTime_Hour_surf
  
  
  ScanTime_Month <- ncvar_get(netcdf, attributes(netcdf$var)$names[24])
  ScanTime_Month_surf <- ScanTime_Month[round(grid_point/49)]
  rm("ScanTime_Month")
  output_matrix[176,23] <- ScanTime_Month_surf
  
  
  flagSurfaceSnowfall <- ncvar_get(netcdf, attributes(netcdf$var)$names[25])
  flagSurfaceSnowfall_surf <- flagSurfaceSnowfall[grid_point]
  rm("flagSurfaceSnowfall")
  output_matrix[176,24] <- flagSurfaceSnowfall_surf
  
  
  binMixedPhaseTop <- ncvar_get(netcdf, attributes(netcdf$var)$names[26])
  binMixedPhaseTop_surf <- binMixedPhaseTop[grid_point]
  rm("binMixedPhaseTop")
  output_matrix[176,25] <- binMixedPhaseTop_surf
  
  
  surfaceSnowfallIndex <- ncvar_get(netcdf, attributes(netcdf$var)$names[27])
  surfaceSnowfallIndex_surf <- surfaceSnowfallIndex[grid_point]
  rm("surfaceSnowfallIndex")
  output_matrix[176,26] <- surfaceSnowfallIndex_surf
  
  
  DSD_phase <- ncvar_get(netcdf, attributes(netcdf$var)$names[28])
  DSD_phase_prof <- vector()
  
  for (i in 1:176) {
    
    DSD_phase_prof[i] <- as.vector(DSD_phase[i,,])[grid_point]
    
  }
  
  rm("DSD_phase")
  output_matrix[,27] <- DSD_phase_prof
  
  
  long_surf <- long[grid_point]
  rm("long")
  output_matrix[176,28] <- long_surf
  
  
  lat_surf <- lat[grid_point]
  rm("lat")
  output_matrix[176,29] <- lat_surf
  
  
  elevation <- ncvar_get(netcdf, attributes(netcdf$var)$names[5])
  elevation_surf <- elevation[grid_point]
  rm("elevation")
  output_matrix[176,30] <- elevation_surf
  
  
  colnames(output_matrix) <- c("airTemperature", "binZeroDeg", "binClutterFreeBottom", "zFactorMeasured_ka", "flagPrecip", 
                               "flagBB", "binBBBottom", "binBBTop", "binBBPeak", "precipRateNearSurface", "precipRate", "zFactorFinal_ka",
                               "paramDSD", "zFactorFinalESurface", "precipRateAve24", "zFactorFinalNearSurface", "precipRateESurface",
                               "phaseNearSurface", "ScanTime_Minute", "ScanTime_Year", "ScanTime_Second", "ScanTime_Hour", 
                               "ScanTime_Month", "flagSurfaceSnowfall", "binMixedPhaseTop", "surfaceSnowfallIndex", "DSD_phase",
                               "long", "lat", "elevation", "zFactorMeasured_ku", "zFactorFinal_ku")
  
  
  write.csv(output_matrix, file = paste0("Directory to save files with profiles of each overpass/GPM_2ADPR_09_Powerline_Vars_KaKu_", substr(files[j],24,39), ".csv"))
  
}


#Generate a timeseries matrix with data from all overpasses

setwd("Directory with files that were processed above")


csv_list <- list.files(getwd())

matrix_ts <- matrix(nrow = 315, ncol = 29)


for (i in 1:315) {
  
  csv_file <- read.csv(csv_list[i])
  
  
  if (i == 290) {
    
    matrix_ts[i,] <- rep(NA, 29)
    
  } else {
    
    
    clutter_free_bin <- csv_file[176,4]
    
    ka_Zm_500m_acfb <- csv_file[clutter_free_bin:(clutter_free_bin - 3),5]
    ka_Zm_500m_acfb <- ifelse(ka_Zm_500m_acfb < 0, NA, ka_Zm_500m_acfb)
    
    ku_Zm_500m_acfb <- csv_file[clutter_free_bin:(clutter_free_bin - 3),32]
    ku_Zm_500m_acfb <- ifelse(ku_Zm_500m_acfb < 0, NA, ku_Zm_500m_acfb)
    
    dfr <- ku_Zm_500m_acfb - ka_Zm_500m_acfb
    ka_Zm_500m_acfb <- ifelse(dfr < -0.5, NA, ka_Zm_500m_acfb)
    ku_Zm_500m_acfb <- ifelse(dfr < -0.5, NA, ku_Zm_500m_acfb)
    
    dfr
    ka_Zm_500m_acfb
    ku_Zm_500m_acfb
    
    
    dfr_cont_flag <- sum(!is.na(ka_Zm_500m_acfb))
    
    ka_cont_flag <- ifelse(sum(is.na(ka_Zm_500m_acfb)) > 1, 0, 1)
    ku_cont_flag <- ifelse(sum(is.na(ku_Zm_500m_acfb)) > 1, 0, 1)
    
    ka_threshold_flag <- ifelse(mean(ka_Zm_500m_acfb, na.rm = T) > 10, 1, 0)
    ku_threshold_flag <- ifelse(sum(ku_Zm_500m_acfb > 8) >= 3, 1, 0)
    
    
    matrix_ts[i,1] <- paste0(csv_file$ScanTime_Year[176], "-", csv_file$ScanTime_Month[176], "-", substr(csv_list[i], 40, 41), " ",
                             csv_file$ScanTime_Hour[176], ":", csv_file$ScanTime_Minute[176], ":", csv_file$ScanTime_Second[176])
    
    matrix_ts[i,2] <- csv_file$lat[176]
    
    matrix_ts[i,3] <- csv_file$long[176]
    
    matrix_ts[i,4] <- csv_file$elevation[176]
    
    matrix_ts[i,5] <- ka_cont_flag
    
    matrix_ts[i,6] <- ka_threshold_flag
    
    matrix_ts[i,7] <- ku_cont_flag
    
    matrix_ts[i,8] <- ku_threshold_flag
    
    matrix_ts[i,9] <- dfr_cont_flag
    
    matrix_ts[i,10] <- ka_Zm_500m_acfb[1]
    
    matrix_ts[i,11] <- mean(ka_Zm_500m_acfb, na.rm = T)
    
    matrix_ts[i,12] <- ku_Zm_500m_acfb[1]
    
    matrix_ts[i,13] <- mean(ku_Zm_500m_acfb, na.rm = T)
    
    matrix_ts[i,14] <- csv_file[clutter_free_bin,13]
    
    matrix_ts[i,15] <- mean(csv_file[clutter_free_bin:(clutter_free_bin - 3),13])
    
    matrix_ts[i,16] <- csv_file[clutter_free_bin,12]
    
    matrix_ts[i,17] <- mean(csv_file[clutter_free_bin:(clutter_free_bin - 3),12])
    
    matrix_ts[i,18] <- csv_file[clutter_free_bin,14]
    
    matrix_ts[i,19] <- mean(csv_file[clutter_free_bin:(clutter_free_bin - 3),14])
    
    matrix_ts[i,20] <- csv_file$zFactorFinalNearSurface[176]
    
    matrix_ts[i,21] <- csv_file$zFactorFinalESurface[176]
    
    matrix_ts[i,22] <- csv_file$precipRateNearSurface[176]
    
    matrix_ts[i,23] <- csv_file$precipRateESurface[176]
    
    matrix_ts[i,24] <- csv_file$precipRateAve24[176]
    
    matrix_ts[i,25] <- csv_file$phaseNearSurface[176]
    
    matrix_ts[i,26] <- clutter_free_bin
    
    matrix_ts[i,27] <- csv_file$binZeroDeg[176]
    
    matrix_ts[i,28] <- csv_file[clutter_free_bin,33]
    
    matrix_ts[i,29] <- mean(csv_file[clutter_free_bin:(clutter_free_bin - 3),33])
    
  }
  
}


colnames(matrix_ts) <- c("datetime", "lat", "long", "elev", "ka_cont_flag", "ka_threshold_flag", "ku_cont_flag", "ku_threshold_flag",
                         "dfr_cont_flag", "nearest_ka_km", "mean_ka_Zm_500m", "nearest_ku_km", "mean_ku_Zm_500m", "nearest_zFactorFinal_ka",
                         "mean_500m_zFactorFinal_ka", "nearest_precipRate", "mean_500m_precipRate", "nearest_paramDSD", "mean_500m_paramDSD",
                         "zFactorFinalNearSurface", "zFactorFinalESurface", "precipRateNearSurface", "precipRateESurface", 
                         "precipRateAve24", "phaseNearSurface", "clutter_free_bin", "binZeroDeg", "nearest_zFactorFinal_ku", 
                         "mean_500m_zFactorFinal_ku")

write.csv(matrix_ts, file = paste0("CSV file to save timeseries matrix"))


#Generate combined file with GPM-DPR data and ground-based data

ground_data <- read.csv("CSV file with ground-based data")[,-c(1)]


complete_matches <- matrix(nrow = 315, ncol = 24)

for (i in 1:315) {
  
  if (i == 290) {
    
    complete_matches[i,1] <- matrix_ts[i,1]
    
    complete_matches[i,2:5] <- rep(NA, 4)
    
    complete_matches[i,6:24] <- as.numeric(ground_data[((which(as.POSIXct(ground_data[,1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT") == round(as.POSIXct(matrix_ts[i,1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), "mins"))) - 360), 2:20])
    
  } else {
    
    complete_matches[i,1] <- matrix_ts[i,1]
    
    complete_matches[i,2] <- matrix_ts[i,11]
    
    complete_matches[i,3] <- matrix_ts[i,13]
    
    complete_matches[i,4] <- matrix_ts[i,15]
    
    complete_matches[i,5] <- matrix_ts[i,29]
    
    complete_matches[i,6:24] <- as.numeric(ground_data[((which(as.POSIXct(ground_data[,1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT") == round(as.POSIXct(matrix_ts[i,1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), "mins"))) - 360), 2:20])
    
  }
  
}

colnames(complete_matches) <- c("timestamp_GMT", "ka_Zm_500m_dBZ", "ku_Zm_500m_dBZ", "ka_zF_500m_dBZ", "ku_zF_500m_dBZ", "air_temp_C", "rel_hum", "wspd_ms", "precip_uncor_mm15", "rain_fraction", "rainfall_mm15", "snowfall_mm15", "hydrometeor_temp_C", "precip_cor_mm15", "precip_type", "a_dry", "a_wet", "part_vel_ms", "part_N0", "part_diam", "Ze_surf_dBz", "Ze_bin_dBZ", "W_surf_ms", "W_bin_ms")


write.csv(complete_matches, file = "CSV file to save combined GPM-DPR and ground-based data")


#Ground clutter analysis

setwd("Directory with GPM-DPR processed data files for each overpass")

gpm_profiles_list <- list.files(getwd())

flag_precip_vec <- vector()

for (i in 1:315) {
  
  if (i == 290) {
    
    flag_precip_vec[i] <- NA
    
  } else {
    
    gpm_profile_matrix <- read.csv(gpm_profiles_list[i])
    
    flag_precip_vec[i] <- gpm_profile_matrix$flagPrecip[176]
    
  }
  
}

length(which(flag_precip_vec != 0))
