#Code to match MRR and GPM-DPR retrievals at different atmospheric levels
#By Andre Bertoncini
#Centre for Hydrology - University of Saskatchewan


library(lubridate)

setwd("set working directory")


#Load MRR data

Ze_files <- list.files(path = getwd(), pattern = "Ze_.*.csv$")

W_files <- list.files(path = getwd(), pattern = "W_.*.csv$")

Time_files <- list.files(path = getwd(), pattern = "Time_.*.csv$")


#Load matrix with data extracted from GPM-DPR files

matrix_ts <- read.csv("CSV file timeseries matrix with GPM-DPR data.csv")[,-c(1)]


#Match MRR profiles with GPM-DPR overpass time

profile_matrix <- matrix(nrow = 31, ncol = 315)

profile_matrix_W <- matrix(nrow = 31, ncol = 315)

for (i in 1:315) {
  
  dates <- as.POSIXct(matrix_ts[i,1], format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
  
  dates <- round(dates, "mins")
  
  k <- which(Ze_files == paste0("Ze_", substr(dates,1,4), substr(dates,6,7), substr(dates,9,10), ".csv"))
  
  if (length(k) == 1) {
    
    Ze <- read.csv(Ze_files[k], header = F)
    
    W <- read.csv(W_files[k], header = F)
    
    time <- as.vector(read.csv(Time_files[k], header = F))
    
    MRR_timestamp <- unname(unlist(time))
    
    MRR_timestamp_gmt <- as.POSIXct(MRR_timestamp, origin = "1970-01-01", tz = "GMT")
    
    Ze_profile <- as.numeric(Ze[which(MRR_timestamp_gmt == dates),])
    
    W_profile <- as.numeric(W[which(MRR_timestamp_gmt == dates),])
    
    profile_matrix[,i] <- rev(Ze_profile)
    
    profile_matrix_W[,i] <- rev(W_profile)
    
  } else {
    
    profile_matrix[,i] <- rep(200, 31)
    
    profile_matrix_W[,i] <- rep(200, 31)
    
  }
  
}


#Filter MRR noise that happens at highest bins during non-precipitating times

profile_matrix <- ifelse(profile_matrix > 100, NA, profile_matrix)

profile_matrix <- ifelse(profile_matrix < 0, NA, profile_matrix)

profile_matrix_W <- ifelse(profile_matrix_W > 100, NA, profile_matrix_W)

profile_matrix_W <- ifelse(profile_matrix_W < -20, NA, profile_matrix_W)

profiles_filtered <- matrix(nrow = 31, ncol = 315)

profiles_filtered_W <- matrix(nrow = 31, ncol = 315)

for (j in 1:315) {
  
  if ((sum(!is.na(profile_matrix[1:11,j])) > 1 & sum(!is.na(profile_matrix[1:11,j])) <= 4) & (sum(is.na(profile_matrix[20:29,j])) == 10)) {
    
    profiles_filtered[1:11,j] <- rep(NA, 11)
    
    profiles_filtered[12:31,j] <- profile_matrix[12:31,j]
    
    profiles_filtered_W[1:11,j] <- rep(NA, 11)
    
    profiles_filtered_W[12:31,j] <- profile_matrix_W[12:31,j]
    
    
  } else {
    
    profiles_filtered[,j] <- profile_matrix[,j]
    
    profiles_filtered_W[,j] <- profile_matrix_W[,j]
    
  }
  
}


#Number of profiles

precipitating_vector <- colSums(!is.na(profiles_filtered))

precipitating_vector <- ifelse(precipitating_vector >= 2, 1, 0)

sum(precipitating_vector)


#Replace profile with NAs during non-precipitating times

for (i in 1:315) {
  
  if (precipitating_vector[i] == 0) {
    
    profiles_filtered[,i] <- rep(NA, 31)
    
    profiles_filtered_W[,i] <- rep(NA, 31)
    
  }
  
}


#Match MRR and GPM-DPR profiles

Ka_Zm <- ifelse(!is.na(matrix_ts[-mrr_na,10]), 1, 0)

Ku_Zm <- ifelse(!is.na(matrix_ts[-mrr_na,12]), 1, 0)


profiles_nona <- ifelse(!is.na(profiles_filtered), 1, 0)

MRR_profiles <- ifelse(colSums(profiles_nona) > 1, 1, 0)


Ka_zFactor <- ifelse(!is.na(matrix_ts[-mrr_na,14]), 1, 0)

Ku_zFactor <- ifelse(!is.na(matrix_ts[-mrr_na,28]), 1, 0)


clutter_free_height <- 22 - (as.numeric(matrix_ts[-mrr_na,26])*0.125)

Ku_matches <- which(Ku_Zm == 1 & MRR_profiles == 1)


#List of files containing data from each GPM-DPR overpass

files_gpm <- list.files("Folder with processed GPM-DPR overpass data")


#Remove data from when MRR had 35-m resolution

profiles_filtered[,1:13] <- matrix(nrow = 31, ncol = 13)

profiles_filtered_W[,1:13] <- matrix(nrow = 31, ncol = 13)


#Retrieve matching retrievals at coincident atmospheric levels

gpm_ka_profile <- list()

gpm_ku_profile <- list()

mrr_profile_list <- list()

mrr_profile_W_list <- list()

gpm_ka_zf_profile <- list()

gpm_ku_zf_profile <- list()

gpm_temp_profile <- list()

precip_rate_profile <- list()


last_index <- c(rep(150,13), rep(111,52), rep (135,250))


for (i in 1:315) {
  
  if (i == 290) {
    
    gpm_ka_profile[[i]] <- NA
    
    gpm_ku_profile[[i]] <- NA
    
    mrr_profile_list[[i]] <- NA
    
    mrr_profile_W_list[[i]] <- NA
    
    gpm_ka_zf_profile[[i]] <- NA
    
    gpm_ku_zf_profile[[i]] <- NA
    
    gpm_temp_profile[[i]] <- NA
    
    precip_rate_profile[[i]] <- NA
    
  } else {
    
    GPM_profiles <- read.csv(paste0("/media/project/abertoncini/02_GPM_MRR/11_GPM_DPR_Data/04_GPM_DPR_Processed_KaKu/", 
                                    files_gpm[i]))
    
    if (GPM_profiles[176,6] == 0) {
      
      Ka_profile <- rep(NA,176)
      Ku_profile <- rep(NA,176)
      Ka_zFactor <- rep(NA,176)
      Ku_zFactor <- rep(NA,176)
      
    } else {
      
      
      Ka_profile <- GPM_profiles[,5]
      Ka_profile <- ifelse(Ka_profile < 0, NA, Ka_profile)
      
      Ku_profile <- GPM_profiles[,32]
      Ku_profile <- ifelse(Ku_profile < 0, NA, Ku_profile)
      
      
      Ka_zFactor <- GPM_profiles[,13]
      Ku_zFactor <- GPM_profiles[,33]
      
    }
    
    
    gpm_temp <- GPM_profiles[,2]
    
    precip_rate <- GPM_profiles[,12]
    
    
    if (i < 66) {
      
      clutter_free_mrr_bin <- 32 - round(((22000 - as.numeric(matrix_ts[i,26])*125)-2141)/200)
      
      
      last_gpm_bin <- vector()
      
      for (j in 0:((as.numeric(matrix_ts[i,26])-135) + 1)) {
        
        last_gpm_bin[j] <- 32 - round((22000 - ((as.numeric(matrix_ts[i,26]) - j)*125)-2141)/200)
        
      }
      
    } else {
      
      
      clutter_free_mrr_bin <- 32 - round(((22000 - as.numeric(matrix_ts[i,26])*125)-2141)/100)
      
      
      last_gpm_bin <- vector()
      
      for (j in 0:((as.numeric(matrix_ts[i,26])-135) + 1)) {
        
        last_gpm_bin[j] <- 32 - round((22000 - ((as.numeric(matrix_ts[i,26]) - j)*125)-2141)/100)
        
      }
      
    }
    
    delta_gpm <- (as.numeric(matrix_ts[i,26]) - (as.numeric(matrix_ts[i,26]) + length(last_gpm_bin)) + 1)
    
    gpm_ka_profile[[i]] <- rev(Ka_profile[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    gpm_ku_profile[[i]] <- rev(Ku_profile[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    gpm_ka_zf_profile[[i]] <- rev(Ka_zFactor[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    gpm_ku_zf_profile[[i]] <- rev(Ku_zFactor[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    gpm_temp_profile[[i]] <- rev(gpm_temp[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    precip_rate_profile[[i]] <- rev(precip_rate[(as.numeric(matrix_ts[i,26]) + delta_gpm):as.numeric(matrix_ts[i,26])])
    
    mrr_profile_list[[i]] <- profiles_filtered[last_gpm_bin,i]
    
    mrr_profile_W_list[[i]] <- profiles_filtered_W[last_gpm_bin,i]
    
    
  }
  
}


gpm_ka_profile_vec <- unlist(gpm_ka_profile)

gpm_ku_profile_vec <- unlist(gpm_ku_profile)

mrr_profile_vec <- unlist(mrr_profile_list)

mrr_profile_W_vec <- unlist(mrr_profile_W_list)

gpm_ka_zf_profile_vec <- unlist(gpm_ka_zf_profile)

gpm_ku_zf_profile_vec <- unlist(gpm_ku_zf_profile)

gpm_temp_vec <- unlist(gpm_temp_profile)

precip_rate_vec <- unlist(precip_rate_profile)


#Generate snowfall rates with retrieved reflectivities

dfr_vec <- gpm_ku_profile_vec - gpm_ka_profile_vec

vel_cm_s <- mrr_profile_W_vec*100

vel_cm_s <- ifelse(vel_cm_s < 0, NA, vel_cm_s)

gpm_temp_vec <- ifelse(gpm_temp_vec > 0, NA, gpm_temp_vec)

N0 <- (5.65*(10^5)*exp((-0.107*gpm_temp_vec)))/100000000

a_coef <- ifelse(gpm_temp_vec < -5, ((3.36*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189), ((5.32*(10000/((N0^(2/3))*(vel_cm_s^(5/3)))))*0.189))

b_exp <- 1.67


#GPM-DPR Snowfall evaluation against MRR estimates

precip_mrr <- ((10^(mrr_profile_vec/10))/a_coef)^(1/b_exp)

precip_ka <- ((10^(gpm_ka_profile_vec/10))/a_coef)^(1/b_exp)

precip_ku <- ((10^(gpm_ku_profile_vec/10))/a_coef)^(1/b_exp)

precip_ka_f <- ((10^(gpm_ka_zf_profile_vec/10))/a_coef)^(1/b_exp)

precip_ku_f <- ((10^(gpm_ku_zf_profile_vec/10))/a_coef)^(1/b_exp)


#Zm-Ka vs. MRR

plot(precip_mrr, precip_ka, ylim = c(0,6.5), xlim = c(0,6.5))

mod_p_ka <- summary(lm(precip_mrr ~ precip_ka))

mod_p_ka

round(sqrt(mod_p_ka$r.squared), 2)

bias = mean(precip_ka - precip_mrr, na.rm = T)
round(bias, 3)

rmse = sqrt(mean((precip_ka - precip_mrr)^2, na.rm = T))
round(rmse, 3)

obs <- precip_mrr
mod <- precip_ka

plot_df <- data.frame(obs, mod)


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = rgb(0,0,1)) +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("Zm-Ka Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()


#Zm-Ku vs. MRR

plot(precip_mrr, precip_ku, ylim = c(0,11), xlim = c(0,11))

mod_p_ku <- summary(lm(precip_mrr ~ precip_ku))

mod_p_ku

round(sqrt(mod_p_ku$r.squared), 2)

bias = mean(precip_ku - precip_mrr, na.rm = T)
round(bias, 3)

rmse = sqrt(mean((precip_ku - precip_mrr)^2, na.rm = T))
round(rmse, 3)


obs <- precip_mrr
mod <- precip_ku

plot_df <- data.frame(obs, mod) 


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = rgb(0,0,1)) +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("Zm-Ku Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()


#Zc-Ka vs. MRR

plot(precip_mrr, precip_ka_f, ylim = c(0,11), xlim = c(0,11))

mod_p_ka_f <- summary(lm(precip_mrr ~ precip_ka_f))

mod_p_ka_f

round(sqrt(mod_p_ka_f$r.squared), 2)

bias = mean(precip_ka_f - precip_mrr, na.rm = T)
round(bias, 3)

rmse = sqrt(mean((precip_ka_f - precip_mrr)^2, na.rm = T))
round(rmse, 3)


obs <- precip_mrr
mod <- precip_ka_f

plot_df <- data.frame(obs, mod) 


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = rgb(1,0,0)) +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("Zc-Ka Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()


#Zc-Ku vs. MRR

plot(precip_mrr, precip_ku_f, ylim = c(0,11), xlim = c(0,11))

mod_p_ku_f <- summary(lm(precip_mrr ~ precip_ku_f))

mod_p_ku_f

round(sqrt(mod_p_ku_f$r.squared), 2)

bias = mean(precip_ku_f - precip_mrr, na.rm = T)
round(bias, 3)

rmse = sqrt(mean((precip_ku_f - precip_mrr)^2, na.rm = T))
round(rmse, 3)


obs <- precip_mrr
mod <- precip_ku_f

plot_df <- data.frame(obs, mod) 


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = rgb(1,0,0)) +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("Zc-Ku Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()


#DPR-P vs. MRR

precip_rate_vec <- ifelse(precip_rate_vec <= 0, NA, precip_rate_vec)

plot(precip_mrr, precip_rate_vec, ylim = c(0,6.5), xlim = c(0,6.5))

mod_p <- summary(lm(precip_mrr ~ precip_rate_vec))

mod_p

sqrt(mod_p$r.squared)

bias = mean(precip_rate_vec - precip_mrr, na.rm = T)
bias

rmse = sqrt(mean((precip_rate_vec - precip_mrr)^2, na.rm = T))
rmse

obs <- precip_mrr
mod <- precip_rate_vec

plot_df <- data.frame(obs, mod) 


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("DPR Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()


#CORRA Evaluation (including matching CORRA with MRR estimates)

gpm_temp_matrix <- matrix(nrow = 176, ncol = 315)

for (i in 1:315) {
  
  if (i == 290) {
    
    gpm_temp_matrix[,i] <- rep(NA, 176)
    
  } else{
    
    GPM_profiles <- read.csv(paste0("Path to files with processed GPM-DPR data at each overpass", files_gpm[i]))
    
    gpm_temp_matrix[,i] <- rev(GPM_profiles[,2])
    
  }
}


kilometric_temp <- matrix(nrow = 88, ncol = 315)

for (j in 0:87) {
  
  kilometric_temp[(j*2 + 2)/2,] <- colMeans(gpm_temp_matrix[(j*2 + 1):(j*2 + 2),])
  
}


index_200m <- 88 - c(79:75,75:71,71:67,67:63,63:59,59:55,55)

temp_200m <- kilometric_temp[index_200m,]

temp_200m <- apply(temp_200m, 2, rev)

Ze_200m <- profiles_filtered

W_200m <- profiles_filtered_W


index_100m <- 88 - c(79:76,76:72,72:68,68,67)

temp_100m <- kilometric_temp[index_100m,]

temp_100m <- apply(temp_100m, 2, rev)


Ze_100m <- matrix(nrow = 16, ncol = 315)

for (j in 0:14) {
  
  Ze_100m[(j*2 + 2)/2,] <- colMeans(profiles_filtered[(j*2 + 1):(j*2 + 2),], na.rm = T)
  
}

Ze_100m[16,] <- profiles_filtered[31,]


W_100m <- matrix(nrow = 16, ncol = 315)

for (j in 0:14) {
  
  W_100m[(j*2 + 2)/2,] <- colMeans(profiles_filtered_W[(j*2 + 1):(j*2 + 2),], na.rm = T)
  
}

W_100m[16,] <- profiles_filtered_W[31,]


Ze_100m <- rbind(matrix(ncol = 315, nrow = 15), Ze_100m)

W_100m <- rbind(matrix(ncol = 315, nrow = 15), W_100m)

Ze_merged <- cbind(Ze_200m[,1:65], Ze_100m[,66:315])

W_merged <- cbind(W_200m[,1:65], W_100m[,66:315])


W_merged <- W_merged*100

W_merged <- ifelse(W_merged < 0, NA, W_merged)


kilometric_temp <- cbind(temp_200m[,1:65], rbind(matrix(ncol = 315, nrow = 15), temp_100m)[,66:315])

kilometric_temp <- ifelse(kilometric_temp > 0, NA, kilometric_temp)


kilometric_N0 <- (5.65*(10^5)*exp((-0.107*kilometric_temp)))/100000000

a_coef_km <- ifelse(kilometric_temp < -5, ((3.36*(10000/((kilometric_N0^(2/3))*(W_merged^(5/3)))))*0.189), ((5.32*(10000/((kilometric_N0^(2/3))*(W_merged^(5/3)))))*0.189))

b_exp <- 1.67

precip_mrr_km <- ((10^(Ze_merged/10))/a_coef_km)^(1/b_exp)


CORRA_input <- as.matrix(read.csv("/media/project/abertoncini/02_GPM_MRR/26_CORRA_5km_Processed/CORRA_5km_precipTotRate_20230905.csv")[,-c(1)])

CORRA_input <- ifelse(CORRA_input <= 0, NA, CORRA_input)

precip_CORRA_200m <- CORRA_input[c(79:75,75:71,71:67,67:63,63:59,59:55,55),1:65]

precip_CORRA_100m <- CORRA_input[c(79:76,76:72,72:68,68,67),66:315]

precip_CORRA_100m <- rbind(precip_CORRA_100m, matrix(ncol = 250, nrow = 15))

precip_CORRA_200m <- apply(precip_CORRA_200m, 2, rev)

precip_CORRA_100m <- apply(precip_CORRA_100m, 2, rev)


precip_CORRA <- cbind(precip_CORRA_200m, precip_CORRA_100m)


#CORRA vs. MRR

precip_corra_vector <- as.vector(precip_CORRA)

precip_mrr_km_vector <- as.vector(precip_mrr_km)

plot(precip_mrr_km_vector, ifelse(precip_corra_vector == 0, NA, precip_corra_vector))


corra_mod <- summary(lm(precip_mrr_km_vector ~ ifelse(precip_corra_vector == 0, NA, precip_corra_vector)))

round(sqrt(corra_mod$r.squared), 2)


bias = mean(ifelse(precip_corra_vector == 0, NA, precip_corra_vector) - precip_mrr_km_vector, na.rm = T)
round(bias,3)

rmse = sqrt(mean((ifelse(precip_corra_vector == 0, NA, precip_corra_vector) - precip_mrr_km_vector)^2, na.rm = T))
round(rmse, 3)


obs <- precip_mrr_km_vector
mod <- ifelse(precip_corra_vector == 0, NA, precip_corra_vector)

plot_df <- data.frame(obs, mod) 


tiff("Path to save plot.tiff", width = 3, height = 2.5, units = "in", res = 1000)

ggplot(plot_df, aes(x = obs, y = mod)) +
  geom_point(alpha = 0.3, size = 0.8) +
  theme_classic() +
  geom_smooth(method = lm, se = TRUE, size = 0.5, color = rgb(0,1,0)) +
  geom_abline(slope = 1, intercept = 0, color = "brown1", linetype = "dashed", size = 0.5) +
  xlab("MRR Snowfall Rate [mm/h]") +
  ylab("CORRA Snowfall Rate [mm/h]") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5)) +
  xlim(0,11) +
  ylim(0,11) +
  theme(text = element_text(size = 10)) +
  theme(axis.line = element_line(size = 0.5))


dev.off()