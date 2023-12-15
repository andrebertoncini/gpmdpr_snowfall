#Code to extract CORRA precipitation rates at Powerline station pixel for each GPM overpass
#By Andre Bertoncini
#Centre for Hydrology - University of Saskatchewan


library(hdf5r)

setwd("Path to CORRA HDF5 files")


files <- list.files()

output_matrix <- matrix(ncol = 315, nrow = 88)


for (j in 1:315) {
  
  hdf5_file <- H5File$new(files[j], mode = "r")
  
  dataset = list.datasets(hdf5_file)
  
  dataset
  
  
  precipTotRate <- readDataSet(hdf5_file[["/KuKaGMI/precipTotRate"]])
  
  long <- readDataSet(hdf5_file[["/KuKaGMI/Longitude"]])
  
  lat <- readDataSet(hdf5_file[["/KuKaGMI/Latitude"]])
  
  
  grid_point <- which.min(sqrt(abs(long - round(-115.198195978996,4))^2 + abs(lat - round(50.8259838210696,5))^2))
  
  
  precipTotRate_prof <- vector()
  
  for (i in 1:88) {
    
    precipTotRate_prof[i] <- (as.vector(precipTotRate[i,,])[grid_point])
    
  }
  
  output_matrix[,j] <- precipTotRate_prof
  
  rm("precipTotRate_prof")
  rm("precipTotRate")
  rm("hdf5_file")
  
}

write.csv(output_matrix, file = "CSV file to save extracted data")

