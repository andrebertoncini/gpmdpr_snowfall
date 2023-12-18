#Code to create a hydrometeor type classification matrix based on Parsivel2 particle diameter and fallspeed.
#By Andre Bertoncini
#Centre for Hydrology - University of Saskatchewan


#Parsivel diameters in mm

diameters <- c(0.062, 0.187, 0.312, 0.437, 0.562, 0.687, 0.812, 0.937, 1.062, 1.187, 1.375, 1.625, 1.875, 2.125, 2.375, 2.750, 3.250, 3.750,
               4.250, 4.750, 5.500, 6.500, 7.500, 8.500, 9.500, 11.000, 13.000, 15.000, 17.000, 19.000, 21.500, 24.500)

#Parsivel velocities in m/s

velocities <- c(0.050, 0.150, 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 0.850, 0.950, 1.100, 1.300, 1.500, 1.700, 1.900, 2.200, 2.600, 3.000,
                3.400, 3.800, 4.400, 5.200, 6.000, 6.800, 7.600, 8.800, 10.400, 12.000, 13.600, 15.200, 17.600, 20.800)

diameters <- diameters[3:32]
velocities <- velocities[3:32]

matrix_output <- matrix(ncol = 30, nrow = 30)

for (j in 1:30) {
  
  for (i in 1:30) {
    
    estimated_v <- c((3.78*(diameters[j]^0.67)), (1.3*(diameters[j]^0.66)), (2.14*((diameters[j]*0.1)^0.20)), (1.07*((diameters[j]*0.1)^0.20)))
    
    estimated_v <- ifelse(estimated_v > 20.800, 20.800, estimated_v)
    
    estimated_d <- c(((velocities[i]/3.78)^(1/0.67)), ((velocities[i]/1.3)^(1/0.66)), ((velocities[i]/2.14)^(1/0.20)), ((velocities[i]/1.07)^(1/0.20)))
    
    estimated_d <- ifelse(estimated_d > 24.500, 24.500, estimated_d)
    
    
    round(estimated_v, 2)
    
    round(estimated_d, 3)
    
    
    if (velocities[i] > estimated_v[1]) {
      
      matrix_output[(31-i),j] <- 1
      
    } else if (velocities[i] < estimated_v[4]) {
      
      matrix_output[(31-i),j] <- 4
      
    } else {
      
      dist <- vector()
      
      
      #Rain
      dist[1] <- sqrt((estimated_v[1] - velocities[i])^2 + (estimated_d[1] - diameters[j])^2)
      
      #Graupel
      dist[2] <- sqrt((estimated_v[2] - velocities[i])^2 + (estimated_d[2] - diameters[j])^2)
      
      #Wet Snow
      dist[3] <- sqrt((estimated_v[3] - velocities[i])^2 + (estimated_d[3] - diameters[j])^2)
      
      #Dry Snow
      dist[4] <- sqrt((estimated_v[4] - velocities[i])^2 + (estimated_d[4] - diameters[j])^2)
      
      
      matrix_output[(31-i),j] <- which.min(dist)
      
    }
    
  }
  
}


write.csv(matrix_output, file = "CSV file to save classification matrix")

