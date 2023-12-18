#Code to process MRR-2, Parsivel2, precipitation gauge, and meteorological data to retrieve near-surface Ze-S Relationships
#By Andre Bertoncini
#Centre for Hydrology - University of Saskatchewan

library(lubridate)
library(ggplot2)


#Process Parsivel data at minute intervals

setwd("Set working directory")


#Load particle classification matrix

matrix_class_input <- as.matrix(read.csv("Load CSV file with Parsivel2 classification matrix")[,c(-1)])

matrix_class <- matrix(ncol = 32, nrow = 32)

matrix_class[1:30,3:32] <- matrix_class_input


#Parsivel list of files

parsivel_list <- list.files(path = getwd())

parsivel_list <- parsivel_list[53:(length(parsivel_list) - 2)]


#Parsivel diameters in mm

diameters <- c(0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625, 1.875, 2.125, 2.375, 2.750, 3.250, 3.750,
               4.250, 4.750, 5.500, 6.500, 7.500, 8.500, 9.500, 11.000, 13.000, 15.000, 17.000, 19.000, 21.500, 24.500)

#Parsivel velocities in m/s

velocities <- c(0.050, 0.150, 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 0.850, 0.950, 1.100, 1.300, 1.500, 1.700, 1.900, 2.200, 2.600, 3.000,
                3.400, 3.800, 4.400, 5.200, 6.000, 6.800, 7.600, 8.800, 10.400, 12.000, 13.600, 15.200, 17.600, 20.800)

#Range of each Parsivel diameter bin in mm

diam_width <- c(0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250, 0.250, 0.5, 0.5, 0.5, 0.5, 0.5,
                1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3)


matrix_v <- matrix(rep(velocities, each = 32), ncol = 32, nrow = 32) #velocity matrix in m/s

matrix_d <- matrix(rep(diameters, 32), ncol = 32, nrow = 32) #diameter matrix in mm

matrix_d_w <- matrix(rep(diam_width, 32), ncol = 32, nrow = 32) #diameter range matrix in mm


#Loop through all 1-minute Parsivel files and calculate Particle Size Distribution parameters.

parsivel_time <- vector()
precip_type <- vector()
part_vel <- vector()
part_diam <- vector()
part_N0 <- vector()
part_a_dry <- vector()
part_a_wet <- vector()
part_Dm <- vector()


for (j in 1:length(parsivel_list)) {
  
  parsivel <- read.table(parsivel_list[j], header = F, sep = ",", dec = ".")
  
  if (ncol(parsivel) == 1041) {
    
    parsivel_time[j] <- as.POSIXct(paste0(parsivel$V1, " ", parsivel$V2), format = "%d.%m.%Y %H:%M:%S", tz = "CST")
    
    matrix_d_v <- matrix(as.numeric(parsivel[,18:1041]), ncol = 32, nrow = 32) #Parsivel matrix (rows = diameters; cols = velocities)
    
    
    sample_area <- (54*(1-(matrix_d/(2*30))))*0.0001 #Parsivel sample area varies with diameter
    
    ND_i <- rowSums(matrix_d_v/(sample_area*60*matrix_v*matrix_d_w), na.rm = T) #m^-3.mm^-1 number concentration
    
    Dliq <- (((6*0.007)/pi)^(1/3))*((diameters*0.1)^(2.2/3))
    
    mD_i <- (Dliq^3)/6
    
    Dm <- sum(ND_i*mD_i*diameters)/sum(ND_i*mD_i)
    
    
    max_ND_index <- which.max(ND_i) #find max ND_i
    
    
    #set ND_is == 0 and ND_is smaller than the diameter equivalent to the maximum ND_i
    
    ND_i_nozero <- ifelse(ND_i[max_ND_index:32] == 0, NA, ND_i[max_ND_index:32])
    
    
    #Requires at least 3 data points to retrieve the N0 intercept
    
    if(sum(!is.na(ND_i_nozero)) >= 3) {
      
      exp_model <- lm(log(ND_i_nozero) ~ diameters[max_ND_index:32]) #compute exponential model
      
      N0 <- unname(exp(exp_model$coefficients[1])*(exp(exp_model$coefficients[2]*(diameters[max_ND_index]))))/100000 #retrieve model intercept
      
      
      vel_m_s <- sum(velocities*colSums(matrix_d_v, na.rm = T))/sum(colSums(matrix_d_v, na.rm = T)) #weighted average of velocities in m/s 
      
      dm_mm <- sum(diameters*rowSums(matrix_d_v, na.rm = T))/sum(rowSums(matrix_d_v, na.rm = T)) #weighted average of diameters in mm
      
      
      vel_cm_s <- vel_m_s*100 #convert velocity to cm/s for calculations below
      
      
      #Compute the a coefficient in Ze = aS^b following Rasmussen et al. (2003). The 0.189 multiplier is to convert Z to Ze in the
      #equation following Smith (1984).
      
      a_dry <- (3.36*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189
      
      a_wet <- (5.32*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189
      
      
      #Determine precip type from Ishizaka et al. (2013): Rain = 1; Graupel = 2; and Rasmussen et al. (1999); Wet snow = 3; Dry snow = 4
      #A classification matrix was generated to reflect the above curves.
      
      precip_type[j] <- unname(matrix_class[(33 - which.min(abs(velocities - vel_m_s))), which.min(abs(diameters - dm_mm))])
      
      
      #Save remaining calculated variables to vector
      
      part_vel[j] <- vel_m_s
      
      part_diam[j] <- dm_mm
      
      part_N0[j] <- N0
      
      part_a_dry[j] <- a_dry
      
      part_a_wet[j] <- a_wet
      
      part_Dm[j] <- Dm
      
      
      #if there are not at least data 3 points, save the variables that do not depend on the exponential model    
      
    } else {
      
      vel_m_s <- sum(velocities*colSums(matrix_d_v, na.rm = T))/sum(colSums(matrix_d_v, na.rm = T))
      
      dm_mm <- sum(diameters*rowSums(matrix_d_v, na.rm = T))/sum(rowSums(matrix_d_v, na.rm = T))
      
      vel_cm_s <- vel_m_s*100
      
      
      precip_type[j] <- unname(matrix_class[(33 - which.min(abs(velocities - vel_m_s))), which.min(abs(diameters - dm_mm))])
      
      part_vel[j] <- vel_m_s
      
      part_diam[j] <- dm_mm
      
      part_N0[j] <- NA
      
      part_a_dry[j] <- NA
      
      part_a_wet[j] <- NA
      
      part_Dm[j] <- NA
      
    }
    
  }
  
}


#Create data frame with outputted vectors

parsivel_output <- data.frame(as.POSIXct(parsivel_time, origin = "1970-01-01", tz = "CST"), precip_type, part_vel, part_diam, part_N0, part_a_dry, part_a_wet)

colnames(parsivel_output) <- c("parsivel_time_CST", "precip_type", "part_vel_ms", "part_diam_mm", "part_N0_m3_mm1", "part_a_dry", "part_a_wet")


#Export data frame to CSV file

write.csv(parsivel_output, file = "CSV file to save Parsivel 1-min processed data")


#Process Parsivel data at 10-sec intervals (SPADE variation)

#Load Parsivel data

setwd("Set working directory for 10-sec Parsivel files")

parsivel_list <- list.files(path = getwd())

header <- read.csv(parsivel_list[1])

parsivel_list <- parsivel_list[2:65]


#Loop through all files and calculate Particle Size Distribution parameters.

for (j in 1:64) {
  
  parsivel <- read.table(parsivel_list[j], header = F, sep = ",", dec = ".", skip = 1)
  
  
  parsivel_time <- vector()
  precip_type <- vector()
  part_vel <- vector()
  part_diam <- vector()
  part_N0 <- vector()
  part_a_dry <- vector()
  part_a_wet <- vector()
  
  
  for (i in 1:nrow(parsivel)) {
    
    parsivel_time[i] <- as.POSIXct(substr(parsivel$V1[i], 1, 19), format = "%Y/%m/%d %H:%M:%S", tz = "GMT")
    
    matrix_d_v <- matrix(as.numeric(parsivel[i,80:1103]), ncol = 32, nrow = 32) #Parsivel matrix (rows = diameters; cols = velocities)
    
    
    sample_area <- (54*(1-(matrix_d/(2*30))))*0.0001 #Parsivel sample area varies with diameter
    
    ND_i <- rowSums(matrix_d_v/(sample_area*10*matrix_v*matrix_d_w), na.rm = T) #m^-3.mm^-1 number concentration
    
    max_ND_index <- which.max(ND_i) #find max ND_i
    
    
    #set ND_is == 0 and ND_is smaller than the diameter equivalent to the maximum ND_i
    
    ND_i_nozero <- ifelse(ND_i[max_ND_index:32] == 0, NA, ND_i[max_ND_index:32])
    
    
    #Requires at least 3 data points to retrieve the N0 intercept
    
    if(sum(!is.na(ND_i_nozero)) >= 3) {
      
      
      exp_model <- lm(log(ND_i_nozero) ~ diameters[max_ND_index:32]) #compute exponential model
      
      N0 <- unname(exp(exp_model$coefficients[1])*(exp(exp_model$coefficients[2]*(diameters[max_ND_index]))))/100000 #retrieve model intercept
      
      
      vel_m_s <- sum(velocities*colSums(matrix_d_v, na.rm = T))/sum(colSums(matrix_d_v, na.rm = T)) #weighted average of velocities in m/s
      
      dm_mm <- sum(diameters*rowSums(matrix_d_v, na.rm = T))/sum(rowSums(matrix_d_v, na.rm = T)) #weighted average of diameters in mm
      
      
      vel_cm_s <- vel_m_s*100 #convert velocity to cm/s for calculations below
      
      
      #Compute the a coefficient in Ze = aS^b following Rasmussen et al. (2003). The 0.189 multiplier is to convert Z to Ze in the
      #equation following Smith (1984).
      
      a_dry <- (3.36*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189
      
      a_wet <- (5.32*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189
      
      
      #Determine precip type from Ishizaka et al. (2013): Rain = 1; Graupel = 2; and Rasmussen et al. (1999); Wet snow = 3; Dry snow = 4
      #A classification matrix was generated to reflect the above curves.
      
      precip_type[i] <- unname(matrix_class[(33 - which.min(abs(velocities - vel_m_s))), which.min(abs(diameters - dm_mm))])
      
      
      #Save remaining calculated variables to vector
      
      part_N0[i] <- N0
      
      part_vel[i] <- vel_m_s
      
      part_diam[i] <- dm_mm
      
      part_a_dry[i] <- a_dry
      
      part_a_wet[i] <- a_wet
      
      
    } else {
      
      precip_type[i] <- NA
      
      part_N0[i] <- NA
      
      part_vel[i] <- NA
      
      part_diam[i] <- NA
      
      part_a_dry[i] <- NA
      
      part_a_wet[i] <- NA
      
    }
    
  }
  
  
  #Create data frame with outputted vectors
  
  parsivel_output <- data.frame(as.POSIXct(parsivel_time, origin = "1970-01-01", tz = "CST"), precip_type, part_vel, part_diam, part_N0, part_a_dry, part_a_wet)
  
  colnames(parsivel_output) <- c("parsivel_time_GMT", "precip_type", "part_vel_ms", "part_diam_mm", "part_N0_m3_mm1", "part_a_dry", "part_a_wet")
  
  
  #Export data frame to CSV file
  
  write.csv(parsivel_output, file = paste0("Directory to save processed 10-sec Parsivel files/Parsivel_SPADE_ProcessedVars_", substr(parsivel_output$parsivel_time_GMT[1], 1, 10), ".csv"))
  
}


#Combine SPADE files into one data.frame

setwd("Set working directory where all the above processed files are")

list_parsivel_proc <- list.files(path = getwd(), pattern = "Parsivel_SPADE_ProcessedVars_.*.csv$")


comb_df <- NULL

for (k in 1:length(list_parsivel_proc)) {
  
  openned_df <- read.csv(list_parsivel_proc[k])
  
  comb_df <- rbind(comb_df, openned_df)
  
}


#Create a Parsivel hourly timestamp 

parsivel_1hr_timestamp <- seq(as.POSIXct(comb_df$parsivel_time_GMT[1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"),
                              as.POSIXct(comb_df$parsivel_time_GMT[length(comb_df$parsivel_time_GMT)], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), 3600) 

time_string <- strftime(as.POSIXct(comb_df$parsivel_time_GMT, format = "%Y-%m-%d %H:%M:%S"), format = "%Y-%m-%d %H:%M")


#Mode function to aggregate particle type

Mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x[!is.na(x)], ux)))]
}


#Aggregate SPADE data to minutes

hourly_matrix <- matrix(nrow = length(parsivel_1hr_timestamp), ncol = 6)

for (l in 1:length(parsivel_1hr_timestamp)) {
  
  hourly_index <- which(substr(parsivel_1hr_timestamp[l], 1, 16) == time_string)
  
  hourly_matrix[l,1] <- Mode(comb_df$precip_type[hourly_index])
  
  hourly_matrix[l,2] <- mean(comb_df$part_vel_ms[hourly_index], na.rm = T)
  
  hourly_matrix[l,3] <- mean(comb_df$part_diam_mm[hourly_index], na.rm = T)
  
  hourly_matrix[l,4] <- mean(comb_df$part_N0_m3_mm1[hourly_index], na.rm = T)
  
  hourly_matrix[l,5] <- mean(comb_df$part_a_dry[hourly_index], na.rm = T)
  
  hourly_matrix[l,6] <- mean(comb_df$part_a_wet[hourly_index], na.rm = T)
  
}  


#Output combined SPADE data

output_SPADE_data <- data.frame(parsivel_1hr_timestamp, hourly_matrix)

colnames(output_SPADE_data) <- c("parsivel_time_GMT", "precip_type", "part_vel_ms", "part_diam_mm", "part_N0_m3_mm1", "part_a_dry", "part_a_wet")

write.csv(output_SPADE_data, file = "CSV file to save combined SPADE data")


#Load precipitation and auxiliary meteorological data

pwl_meteo <- read.table("TXT file with gauge precipitation data and auxiliary meteo data", header = T, sep = ",", dec = ".", skip = 1)

pwl_meteo <- pwl_meteo[32:330367,]


pwl_ta <- as.numeric(pwl_meteo$Air.temperature)
pwl_ta <- ifelse(pwl_ta == -9999, NA, pwl_ta)
plot(pwl_ta)

pwl_rh <- as.numeric(pwl_meteo$Relative.Humidity)
pwl_rh <- ifelse(pwl_rh == -9999, NA, pwl_rh)
plot(pwl_rh)

pwl_wsp <- as.numeric(pwl_meteo$Wind.Speed.B)
pwl_wsp <- ifelse(pwl_wsp == -9999, NA, pwl_wsp)
plot(pwl_wsp)

pwl_precip <- as.numeric(pwl_meteo$Incremental.All.Precipitation.Pluvio)
pwl_precip <- ifelse(pwl_precip == -9999, NA, pwl_precip)
plot(pwl_precip)


#Correct for snowfall undercatch

#Precipitation partitioning (Harder and Pomeroy, 2013)

P <- pwl_precip
T_C <- pwl_ta
RHl <- pwl_rh
T_K <- T_C + 273.15


M_w <- 0.01801528 #kg/mol
R <- 8.31441 #J/mol.K

D = 2.06e-5*(((T_K)/273.15)^1.75)

lambda_t = (0.000063*(T_K)) + 0.00673

L = ifelse(T_C < 0, 1000*(2834.1 - (0.29*T_C) - (0.004*(T_C^2))), 1000*(2501 - 2.36*T_C))

e = RHl/100*(0.611*(exp((17.3*T_C)/(237.3+T_C))))

rho_Ta = ((M_w*e)/(R*T_K))*1e+3

#iterative process

T_i = 271.5 #initial guess

for (i in 100) {
  
  e_ti = RHl/100*(0.611*(exp((17.3*(T_i - 273.15))/(237.3+(T_i - 273.15)))))
  
  es_ti = e_ti*100/RHl
  
  rho_Ti = ((M_w*es_ti)/(R*T_i))*1e+3
  
  T_i = T_K + ((D/lambda_t)*L)*(rho_Ta - rho_Ti)
  
}

#Fraction of precipitation as rainfall

T_fr <- T_i-273.15

b <- 2.630006
c <- 0.09336

Fr = 1/(1+b*(c^T_fr))

Rainfall = P*Fr
Snowfall = P*(1-Fr)


Ws <- pwl_wsp

CE = 1.18*exp(-0.18*Ws)

CE <- ifelse(CE > 1, 1, CE)

Snowfall_Cor <- Snowfall/CE

Precip_Cor <- Snowfall_Cor + Rainfall

plot(Precip_Cor)


precip_timestamp_cst_15 <- seq(as.POSIXct(pwl_meteo$TIMESTAMP[1], format = "%Y-%m-%d %H:%M:%S", tz = "CST"),
                               as.POSIXct(pwl_meteo$TIMESTAMP[length(pwl_meteo$TIMESTAMP)], format = "%Y-%m-%d %H:%M:%S", tz = "CST"), 900) 


meteo_output <- data.frame(as.character(precip_timestamp_cst_15), pwl_ta, pwl_rh, pwl_wsp, pwl_precip, Fr, Rainfall, Snowfall, T_fr, Precip_Cor)

colnames(meteo_output) <- c("time_CST", "ta", "rh", "wspd", "precip_uncor", "rain_fraction", "rainfall", "snowfall", "hydro_temp", "precip_cor")

write.csv(meteo_output, file = "CSV files to save processed meteo outputs")



#Create hourly timestamp for precipitation data

precip_timestamp_cst <- seq(as.POSIXct(pwl_meteo$TIMESTAMP[1], format = "%Y-%m-%d %H:%M:%S", tz = "CST"),
                            as.POSIXct(pwl_meteo$TIMESTAMP[length(pwl_meteo$TIMESTAMP)], format = "%Y-%m-%d %H:%M:%S", tz = "CST"), 3600) 


#Reload Parsivel data from SPADE

Parsivel_SPADE <- read.csv("CSV file with hourly SPADE data")

parsivel_spade_timestamp_CST <- seq(as.POSIXct("2019-04-23 18:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "CST"),
                                    as.POSIXct("2019-06-26 11:10:00", format = "%Y-%m-%d %H:%M:%S", tz = "CST"), 3600)

Parsivel_SPADE_CST <- data.frame(parsivel_spade_timestamp_CST, Parsivel_SPADE[,2:8])


#Reload my Parsivel data

My_Parsivel <- read.csv("CSV file with 1-min processed parsivel data")


parsivel_data_full <- My_Parsivel[!is.na(My_Parsivel$precip_type),]

parsivel_timestamp_cst <- as.POSIXct(parsivel_data_full$parsivel_time_CST, format = "%Y-%m-%d %H:%M:%S", tz = "CST")


#Aggregate precipitation to hourly

len = length(precip_timestamp_cst)
start = 1
end = 4

precip_mm_hr <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  precip_sub <- Precip_Cor[start:end]
  precip_sum <- sum(precip_sub, na.rm = T)
  precip_mm_hr[i] <- precip_sum
  
  start = start+4
  end = end+4
}


#Aggregate hydrometeor temperature to hourly

len = length(precip_timestamp_cst)
start = 1
end = 4

T_h_hr <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  T_h_sub <- T_fr[start:end]
  T_h_mean <- mean(T_h_sub, na.rm = T)
  T_h_hr[i] <- T_h_mean
  
  start = start+4
  end = end+4
}


#Aggregate air temperature to hourly

len = length(precip_timestamp_cst)
start = 1
end = 4

T_a_hr <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  T_a_sub <- T_C[start:end]
  T_a_mean <- mean(T_a_sub, na.rm = T)
  T_a_hr[i] <- T_a_mean
  
  start = start+4
  end = end+4
}


#Match precipitation and My Parsivel data

precip_df <- data.frame(precip_timestamp_cst, precip_mm_hr, T_a_hr, T_h_hr)

colnames(precip_df) <- c("time_CST", "precip_cor", "ta", "hydro_temp")


parsivel_minute_cst <- seq(as.POSIXct(precip_df$time_CST[1], format = "%Y-%m-%d %H:%M:%S", tz = "CST"),
                           as.POSIXct(precip_df$time_CST[length(precip_df$time_CST)], format = "%Y-%m-%d %H:%M:%S", tz = "CST"), 60)


parsivel_df_min <- data.frame(parsivel_minute_cst, rep(NA, length(parsivel_minute_cst)), rep(NA, length(parsivel_minute_cst)), rep(NA, length(parsivel_minute_cst)), rep(NA, length(parsivel_minute_cst)), rep(NA, length(parsivel_minute_cst)), rep(NA, length(parsivel_minute_cst)))

colnames(parsivel_df_min) <- c("parsivel_minute_cst", "precip_type", "part_a_dry", "part_a_wet", "part_vel", "part_N0", "part_diam")

parsivel_df_min$precip_type[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$precip_type

parsivel_df_min$part_a_dry[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$part_a_dry

parsivel_df_min$part_a_wet[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$part_a_wet

parsivel_df_min$part_vel[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$part_vel_ms

parsivel_df_min$part_N0[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$part_N0_m3_mm1

parsivel_df_min$part_diam[which(parsivel_minute_cst %in% parsivel_timestamp_cst)] <- parsivel_data_full$part_diam_mm


#Aggregate My Parsivel data to hourly

Mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x[!is.na(x)], ux)))]
}


len = length(precip_df$time_CST)
start = 1
end = 60

precip_type_15min <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  type_sub <- parsivel_df_min$precip_type[start:end]
  type_mode <- Mode(type_sub)
  precip_type_15min[i] <- type_mode
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

a_dry_h <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  a_dry_sub <- parsivel_df_min$part_a_dry[start:end]
  a_dry_mean <- mean(a_dry_sub, na.rm = T)
  a_dry_h[i] <- a_dry_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

a_wet_h <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  a_wet_sub <- parsivel_df_min$part_a_wet[start:end]
  a_wet_mean <- mean(a_wet_sub, na.rm = T)
  a_wet_h[i] <- a_wet_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

vel_h <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  vel_sub <- parsivel_df_min$part_vel[start:end]
  vel_mean <- mean(vel_sub, na.rm = T)
  vel_h[i] <- vel_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

diam_h <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  diam_sub <- parsivel_df_min$part_diam[start:end]
  diam_mean <- mean(diam_sub, na.rm = T)
  diam_h[i] <- diam_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

N0_h <- vector(mode = "numeric", length = length(precip_df$time_CST))

for (i in 1:len) {
  N0_sub <- parsivel_df_min$part_N0[start:end]
  N0_mean <- mean(N0_sub, na.rm = T)
  N0_h[i] <- N0_mean
  
  start = start+60
  end = end+60
}


precip_parsivel_df <- data.frame(precip_df, precip_type_15min, a_dry_h, a_wet_h, vel_h, diam_h, N0_h)


#Load MRR data

setwd("Set working directory with MRR-2 CSV files")


Ze_files <- list.files(path = getwd(), pattern = "Ze_.*.csv$")

Time_files <- list.files(path = getwd(), pattern = "Time_.*.csv$")

W_files <- list.files(path = getwd(), pattern = "W_.*.csv$")

Ze_list_surf <- list()
Ze_list_bin <- list()
W_list_surf <- list()
W_list_bin <- list()
Time_list <- list()

bin_number <- c(rep(30,51), rep(7,130), rep(15,946))


#Retrieve Ze and time from MRR-2 data for first return

for (k in 1:length(Ze_files)) {
  
  Ze <- read.csv(Ze_files[k], header = F)
  W <- read.csv(W_files[k], header = F)
  
  Ze <- ifelse(as.matrix(Ze) == -9999, NA, as.matrix(Ze))
  W <- ifelse(as.matrix(W) == -9999, NA, as.matrix(W))
  
  Ze_list_surf[[k]] <- Ze[,3]
  Ze_list_bin[[k]] <- Ze[,bin_number[k]]
  
  W_list_surf[[k]] <- W[,3]
  W_list_bin[[k]] <- W[,bin_number[k]]
  
  
  Time_list[[k]] <- as.vector(read.csv(Time_files[k], header = F))
  
}


Ze_timeseries_surf <- unlist(Ze_list_surf)
Ze_timeseries_bin <- unlist(Ze_list_bin)

W_timeseries_surf <- unlist(W_list_surf)
W_timeseries_bin <- unlist(W_list_bin)

MRR_timestamp <- unname(unlist(Time_list))

MRR_timestamp_gmt <- as.POSIXct(MRR_timestamp, origin = "1970-01-01", tz = "GMT")


#Match MRR and Parsivel

mrr_minute_gmt <- seq(as.POSIXct(precip_df$time_CST[1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"),
                      as.POSIXct(precip_df$time_CST[length(precip_df$time_CST)], format = "%Y-%m-%d %H:%M:%S", tz = "GMT"), 60)

MRR_minute_df <- data.frame(mrr_minute_gmt, rep(NA, length(mrr_minute_gmt)), rep(NA, length(mrr_minute_gmt)), rep(NA, length(mrr_minute_gmt)), rep(NA, length(mrr_minute_gmt))) 

colnames(MRR_minute_df) <- c("mrr_minute_gmt", "MRR_Ze_surf", "MRR_Ze_bin", "MRR_W_surf", "MRR_W_bin")

MRR_minute_df$MRR_Ze_surf[which(mrr_minute_gmt %in% MRR_timestamp_gmt)] <- Ze_timeseries_surf
MRR_minute_df$MRR_Ze_bin[which(mrr_minute_gmt %in% MRR_timestamp_gmt)] <- Ze_timeseries_bin
MRR_minute_df$MRR_W_surf[which(mrr_minute_gmt %in% MRR_timestamp_gmt)] <- W_timeseries_surf
MRR_minute_df$MRR_W_bin[which(mrr_minute_gmt %in% MRR_timestamp_gmt)] <- W_timeseries_bin

MRR_minute_df$MRR_Ze_surf_cst <- c(MRR_minute_df$MRR_Ze_surf[361:length(MRR_minute_df$MRR_Ze_surf)], rep(NA, 360))
MRR_minute_df$MRR_Ze_bin_cst <- c(MRR_minute_df$MRR_Ze_bin[361:length(MRR_minute_df$MRR_Ze_bin)], rep(NA, 360))
MRR_minute_df$MRR_W_surf_cst <- c(MRR_minute_df$MRR_W_surf[361:length(MRR_minute_df$MRR_W_surf)], rep(NA, 360))
MRR_minute_df$MRR_W_bin_cst <- c(MRR_minute_df$MRR_W_bin[361:length(MRR_minute_df$MRR_W_bin)], rep(NA, 360))


#Aggregate MRR to hourly

len = length(precip_df$time_CST)
start = 1
end = 60

ze_hourly_surf <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  ze_sub <- MRR_minute_df$MRR_Ze_surf_cst[start:end]
  ze_mean <- mean(ze_sub, na.rm = T)
  ze_hourly_surf[i] <- ze_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

ze_hourly_bin <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  ze_sub <- MRR_minute_df$MRR_Ze_bin_cst[start:end]
  ze_mean <- mean(ze_sub, na.rm = T)
  ze_hourly_bin[i] <- ze_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

W_hourly_surf <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  ze_sub <- MRR_minute_df$MRR_W_surf_cst[start:end]
  ze_mean <- mean(ze_sub, na.rm = T)
  W_hourly_surf[i] <- ze_mean
  
  start = start+60
  end = end+60
}


len = length(precip_df$time_CST)
start = 1
end = 60

W_hourly_bin <- vector(mode = "numeric", length = len)

for (i in 1:len) {
  ze_sub <- MRR_minute_df$MRR_W_bin_cst[start:end]
  ze_mean <- mean(ze_sub, na.rm = T)
  W_hourly_bin[i] <- ze_mean
  
  start = start+60
  end = end+60
}


precip_parsivel_mrr_df <- data.frame(precip_parsivel_df, ze_hourly_surf, ze_hourly_bin, W_hourly_surf, W_hourly_bin)


#Add SPADE parsivel data full database

precip_parsivel_mrr_df$precip_type[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$precip_type

precip_parsivel_mrr_df$a_dry_h[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$part_a_dry

precip_parsivel_mrr_df$a_wet_h[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$part_a_wet

precip_parsivel_mrr_df$vel_h[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$part_vel_ms

precip_parsivel_mrr_df$N0_h[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$part_N0_m3_mm1

precip_parsivel_mrr_df$diam_h[which(precip_parsivel_mrr_df$time_CST %in% Parsivel_SPADE_CST$parsivel_spade_timestamp_CST)] <- Parsivel_SPADE_CST$part_diam_mm


precip_parsivel_mrr_df <- precip_parsivel_mrr_df[,-c(15)]


#Export data

colnames(precip_parsivel_mrr_df) <- c("timestamp_CST", "precip_cor_mm", "air_temp_C", "hydrometeor_temp_C", "precip_type", "a_dry", "a_wet", "part_vel_ms", "part_diam_mm", "part_N0_m3_mm1", "Ze_surf_dBz", "Ze_bin_dBZ", "W_surf_ms", "W_bin_ms")

write.csv(precip_parsivel_mrr_df, file = "CSV file to save the complete processed file with gauge, meteo, MRR-2, and Parsivel2 data")


#Subset complete database to solid precipitation

precip_parsivel_mrr_df <- read.csv("Reload CSV file to save the complete processed file with gauge, meteo, MRR-2, and Parsivel2 data")

solid_precip_df <- precip_parsivel_mrr_df[which(precip_parsivel_mrr_df$precip_type != 1),]


#Real-time a coefficient for theoretical method

a_precip_type <- ifelse(solid_precip_df$precip_type == 4, solid_precip_df$a_dry, solid_precip_df$a_wet)


#Snowfall evaluation using Parsivel2 observed N0

obs_precip <- solid_precip_df$precip_cor_mm

Ze_mm6mm3 <- 10^(solid_precip_df$Ze_surf_dBz/10)

Ze_mm6mm3 <- Ze_mm6mm3[which(obs_precip != 0)]

obs_precip <- obs_precip[which(obs_precip != 0)]

Ze_mm6mm3 <- 10^(solid_precip_df$Ze_surf_dBz/10)

S_est <- (Ze_mm6mm3/a_precip_type)^(1/1.67)


plot(solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], ylim = c(0,6.5), xlim = c(0,6.5))

round(cor(solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]), 2)


bias <- mean(S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0] - solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0])
round(bias,3)

rmse <- sqrt(mean((S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0] - solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0])^2))
round(rmse,3)


obs <- solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]

mod <- S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]

plot_df <- data.frame(obs, mod) 


tiff("TIFF file to save scatterplot", width = 6, height = 5, units = "in", res = 600)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.8, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed") +
  xlab("Observed Snowfall Rate [mm/h]") +
  ylab("Estimated Snowfall Rate (N0_obs) [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlim(0,6.5) +
  ylim(0,6.5) +
  theme(text = element_text(size = 16))

dev.off()


#Snowfall evaluation using temperature estimated N0

N0_mod <- (5.65*(10^5)*exp((-0.107*solid_precip_df$air_temp_C)))/100000000

a_mod <- ifelse(solid_precip_df$air_temp_C < -5, 3.36e+4/((N0_mod^(2/3))*((solid_precip_df$part_vel_ms*100)^(5/3))), 5.32e+4/((N0_mod^(2/3))*((solid_precip_df$part_vel_ms*100)^(5/3))))*0.189


obs_precip <- solid_precip_df$precip_cor_mm

Ze_mm6mm3 <- 10^(solid_precip_df$Ze_surf_dBz/10)

Ze_mm6mm3 <- Ze_mm6mm3[which(obs_precip != 0)]

obs_precip <- obs_precip[which(obs_precip != 0)]

Ze_mm6mm3 <- 10^(solid_precip_df$Ze_surf_dBz/10)

S_est <- (Ze_mm6mm3/a_mod)^(1/1.67)


plot(solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], ylim = c(0,6.5), xlim = c(0,6.5))

round(cor(solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0], S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]), 2)


bias <- mean(S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0] - solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0])
round(bias,3)

rmse <- sqrt(mean((S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0] - solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0])^2))
round(rmse,3)


obs <- solid_precip_df$precip_cor_mm[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]

mod <- S_est[!is.na(S_est) & solid_precip_df$precip_cor_mm != 0]

plot_df <- data.frame(obs, mod) 


tiff("TIFF file to save scatterplot", width = 6, height = 5, units = "in", res = 600)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.8, color = "purple") +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed") +
  xlab("Observed Snowfall Rate [mm/h]") +
  ylab("Estimated Snowfall Rate (N0_est) [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  xlim(0,6.5) +
  ylim(0,6.5) +
  theme(text = element_text(size = 16))

dev.off()

