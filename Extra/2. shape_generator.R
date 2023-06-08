# load required packages:
library(terra)
library(geodata)
# import a world countries map:
countries <- world(resolution = 5, path = "maps")  # you may choose a smaller (more detailed) resolution for the polygon borders, and a different folder path to save the imported map
head(countries)
# import a table with country codes and continents:
cntry_codes <- country_codes()
head(cntry_codes)
# add this table to the countries map attributes:
head(cntry_codes[ , 1:4])

countries <- merge(countries, cntry_codes, by.x = "GID_0", by.y = "ISO3", all.x = TRUE)

# plot the countries map coloured according to "continent":
plot(countries, "continent", lwd = 0.2, main = "Countries by continent")

# dissolve (aggregate) countries into a continents map:
continents <- aggregate(countries, by = "continent")
values(continents)
plot(continents, "continent", lwd = 0.2)

# note that each continent (colour) is a multi-part polygon including mainland and islands - see also:
plot(continents[1, ])
terra::writeVector(continents[1, ], filename="TESTO3", overwrite=TRUE)


# disaggregate continent polygons, to then separate islands and mainlands:
continents <- disagg(continents)
# get a map of just the continent mainlands (largest polygons):
unique(continents$continent)
largest <- (order(expanse(continents), decreasing = TRUE))[1:length(unique(continents$continent))]
mainlands <- continents[largest, ]
plot(mainlands[1,], "continent", lwd = 0.2, main = "Continent mainlands")



terra::writeVector(mainlands[1, ], filename="TESTO2", overwrite=TRUE)

