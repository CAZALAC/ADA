clipped_virtualstations <- raster::intersect(states,virtualstations)
states <- get_country_shape(country="BWA")
sf_objet <- sf::st_as_sf(geodata::gadm(country=country, level=0, path=tempdir()))
#legacy
sf_objet <- as(sf_objet, "Spatial")
return(sf_objet)

ISO.codes=read.csv("CountryISOCodes.csv",sep=";")
Afr.country.list=as.character(ISO.codes$ThreeLetter)
virtualstations <- raster::shapefile("./Shape/chirps_25_lite_2.shp")

for (country2 in Afr.country.list) {
  print(paste("Working... ",country2))
  states <- get_country_shape(country=country2)
  clipped_virtualstations <- raster::intersect(states,virtualstations)
  shapefile(clipped_virtualstations,paste("./Shape/",country2,"_chirps_25.shp",sep=""), overwrite=TRUE)
}
