#cargamos el archivo con metadata de las estaciones
#metadata <- read.csv("BWA/BaseDatosEstaciones.csv", header = TRUE, sep = ",", dec = ".")
#cargamos la data de las estaciones
data <- read.csv("BWA/BaseDatosRegistros.csv", header = TRUE, sep = ",", dec = ".")
data[3:14] <- data[3:14]*0.7
#guardamos el archivo data en un archivo csv
write.csv(data, file = "BWA/BaseDatosRegistros.csv", row.names = FALSE)