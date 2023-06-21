library(shiny);library(leaflet);library(shinybusy);library(shinyjs);library(raster);library(rgdal);library(countrycode);library(rts);require(ncdf4);library(HelpersMG)
library(plyr);library(latticeExtra);library(reshape2);library(gdata);library(corrplot);library(sqldf)
library(zoo);library(Kendall);library(zyp);library(car);library(gtools);
library(rgeos);library(lmom);library(lmomRFA);library(sp);library(rrcov);library(nsRFA);library(ModelMap)
library(stringr);library(rasterVis);library(hydroGOF);library(randomForest);library(progress);
library(gtools);library(here);library(chron);library(lattice);library(RColorBrewer);
library(sf);library(circular);library(reshape);library(deldir);library(shinyalert);library(geodata)

workdir <- here()
setwd(workdir)

source('DroughAtlasFunctions.R')

country_list <- country_list()
datadownload <- c("CRU","CHIRPS")
resolution <- c("0.25 degrees","0.05 degrees")
# ui object
ui <- fluidPage(
  titlePanel(p("African Drought Atlas", style = "color:#3474A7")),
  
  tabsetPanel(
    tabPanel("Run", 

             sidebarLayout(
               sidebarPanel(
                 fileInput(
                   inputId = "filemap",
                   label = "Upload map. Choose shapefile",
                   multiple = TRUE,
                   accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj")
                 ),
                 tags$h5("or"),
                 selectInput("icountry_list", "Select a Country", country_list),
                 selectInput(inputId = "ddata", "Data Source", datadownload ),
                 conditionalPanel(condition="input.ddata == 'CHIRPS'", 
                                  selectInput(inputId = "resol", "Resolution", resolution ),
                 ),
                 actionButton(inputId = "Runs1", label = "Run Step 1"),
                 actionButton(inputId = "Runs2", label = "Run Steps 2-7"),
                 actionButton(inputId = "RunsA", label = "Run All"),
                 useShinyjs(),
                 actionButton("refresh", "Refresh"),

                 
                 
               ),
               
               mainPanel(
                 leafletOutput("map"),
                 textOutput("stations")
                 
               )
             )           
             
    ),
    tabPanel("Visualise Results", 
             
             sidebarLayout(
               sidebarPanel(
                 fileInput(
                   inputId = "filemap2",
                   label = "Upload map. Choose shapefile",
                   multiple = TRUE,
                   accept = c(".shp", ".dbf", ".sbn", ".sbx", ".shx", ".prj")
                 ),
                 tags$h5("or"),
                 selectInput("icountry_list2", "Select a Country", country_list),
                 actionButton(inputId = "RunsV", label = "Load Data")
                 
               ),
               
               mainPanel(
                 leafletOutput("map2")
                 )
             ) ),
    tabPanel("Config", 
             
             sidebarLayout(
               sidebarPanel(
                 selectInput("clipping", "Clipping Method", c("Rectangle","Shape")),
                 p("The Shape option only occupies stations that are within the area of interest. The rectangle option uses the extend to generate a rectangle with 5% margin."),
                 numericInput("maxstations", "Max. Stations:", 1000, min = 1, max = 10000),
                 p("Max number of virtual stations to generate."),
                 
                 
              ),
               
               mainPanel(
                 #plotOutput("mapPlot"),
                 #plotOutput("mapPlot2"),
                 #leafletOutput("map2")
                 
               )
             ) ),
    

    
    ),

)

# server()
server <- function(input, output) {

  #Get and print the numbers of stations 
  output$stations <- renderText({
    
    # Name or shape detection
    if(length(input$filemap) > 1 ){
      req(input$filemap)
      shpdf <- input$filemap
      tempdirname <- dirname(shpdf$datapath[1])
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      countryiso = shpdf$name[grep(pattern = "*.shp$", shpdf$name)] 
      countryiso = tools::file_path_sans_ext(countryiso)
      
    }else{
      countryiso = countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c')
    }
    
    
  if(file.exists(paste0(countryiso,"/BaseDatosEstaciones.csv", sep=""))){
    try(stations <- read.csv(paste0(countryiso,"/BaseDatosEstaciones.csv", sep="")))
    return(paste("Stations:",length(unique(stations$id_station)), " stations detected!")) 
  } else {
    return("Stations: No data")
  }

  })
  
  observeEvent(input$refresh, {
    refresh()
  })
  #Get country shape
  boundariescountry <- reactive({
    boundariescountry <- get_country_shape(country=input$icountry_list) #raster::getData('GADM', country=countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c'), level=0)
    boundariescountry
  })
  #hard codding, arreglar
  boundariescountry2 <- reactive({
    boundariescountry2 <- get_country_shape(country=input$icountry_list2) #raster::getData('GADM', country=countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c'), level=0)
    boundariescountry2
  })

  #The Default map
  output$map <- renderLeaflet({
    if(length(input$filemap) > 1 )
    {
      req(input$filemap)
      
      # shpdf is a data.frame with the name, size, type and
      # datapath of the uploaded files
      shpdf <- input$filemap
      
      # The files are uploaded with names
      # 0.dbf, 1.prj, 2.shp, 3.xml, 4.shx
      # (path/names are in column datapath)
      # We need to rename the files with the actual names:
      # fe_2007_39_county.dbf, etc.
      # (these are in column name)
      
      # Name of the temporary directory where files are uploaded
      tempdirname <- dirname(shpdf$datapath[1])
      
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      
      # Now we read the shapefile with readOGR() of rgdal package
      # passing the name of the file with .shp extension.
      
      # We use the function grep() to search the pattern "*.shp$"
      # within each element of the character vector shpdf$name.
      # grep(pattern="*.shp$", shpdf$name)
      # ($ at the end denote files that finish with .shp,
      # not only that contain .shp)

      map <- raster::shapefile(paste(tempdirname,
                              shpdf$name[grep(pattern = "*.shp$", shpdf$name)],
                              sep = "/"))
      try(map <- spTransform(map, "+proj=longlat +datum=WGS84"), silent = TRUE)
      
      leaflet() %>%addTiles() %>%
        addPolygons(data= map, color = "#444444", weight = 1, smoothFactor = 0.5,
                    opacity = 0.5, fillOpacity = 0,
                    highlightOptions = highlightOptions(color = "white", weight = 2,
                                                        bringToFront = TRUE), group = "Boundaries")  %>%
  addLayersControl(
    overlayGroups = c("Boundaries"),
    options = layersControlOptions(collapsed = FALSE)
  )
      
    } else {
      
      states <- boundariescountry()

      
      m <- leaflet() %>%addTiles() %>%
        addPolygons(data= states, color = "#444444", weight = 1, smoothFactor = 0.5,
                    opacity = 0.5, fillOpacity = 0,
                    highlightOptions = highlightOptions(color = "white", weight = 2,
                                                        bringToFront = TRUE), group = "Boundaries") %>%
        
        
        addLayersControl(
          overlayGroups = c("Boundaries"),
          options = layersControlOptions(collapsed = FALSE)
        )
      
      m
                                                
    }
  })
  
  
  # Visualization map
  output$map2 <- renderLeaflet({
    if(length(input$filemap2) > 1 )
    {
      
      req(input$filemap2)
      
      # shpdf is a data.frame with the name, size, type and
      # datapath of the uploaded files
      shpdf <- input$filemap2
      
      # The files are uploaded with names
      # 0.dbf, 1.prj, 2.shp, 3.xml, 4.shx
      # (path/names are in column datapath)
      # We need to rename the files with the actual names:
      # fe_2007_39_county.dbf, etc.
      # (these are in column name)
      
      # Name of the temporary directory where files are uploaded
      tempdirname <- dirname(shpdf$datapath[1])
      
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      
      # Now we read the shapefile with readOGR() of rgdal package
      # passing the name of the file with .shp extension.
      
      # We use the function grep() to search the pattern "*.shp$"
      # within each element of the character vector shpdf$name.
      # grep(pattern="*.shp$", shpdf$name)
      # ($ at the end denote files that finish with .shp,
      # not only that contain .shp)
      #map <- raster::shapefile("AfricaDA.shp")
      
      
      
      map <- raster::shapefile(paste(tempdirname,
                                     shpdf$name[grep(pattern = "*.shp$", shpdf$name)],
                                     sep = "/"))
      try(map <- spTransform(map, "+proj=longlat +datum=WGS84"), silent = TRUE)
      
      # map <- spTransform(map, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
      #raster::projection(dfgff) <- 
      #dfgff
      #map <- sf::st_read()
      #transformamos por si viene en utm y dejar en lat lon
      #map <- st_transform(map, "+init=epsg:4326")
      
      
      leaflet() %>%addTiles() %>%
        addPolygons(data= map, color = "#444444", weight = 1, smoothFactor = 0.5,
                    opacity = 0.5, fillOpacity = 0,
                    highlightOptions = highlightOptions(color = "white", weight = 2,
                                                        bringToFront = TRUE), group = "Boundaries")  %>%
        addLayersControl(
          overlayGroups = c("Boundaries"),
          options = layersControlOptions(collapsed = FALSE)
        )
      
    } else {
      
      states <- boundariescountry2()
      
      
      m <- leaflet() %>%addTiles() %>%
        addPolygons(data= states, color = "#444444", weight = 1, smoothFactor = 0.5,
                    opacity = 0.5, fillOpacity = 0,
                    highlightOptions = highlightOptions(color = "white", weight = 2,
                                                        bringToFront = TRUE), group = "Boundaries") %>%
        
        
        addLayersControl(
          overlayGroups = c("Boundaries"),
          options = layersControlOptions(collapsed = FALSE)
        )
      
      m
      
    }
  })
  
  
  
  Output_f2 = reactiveVal()
  Output_f3 = reactiveVal()
  Output_f4 = reactiveVal()
  Output_f5 = reactiveVal()
  
  
  observeEvent(input$Runs1, {
    
    if(length(input$filemap) > 1 )
      
    {
      req(input$filemap)
      shpdf <- input$filemap
      tempdirname <- dirname(shpdf$datapath[1])
      
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      
      countryiso = shpdf$name[grep(pattern = "*.shp$", shpdf$name)] 
      countryiso = tools::file_path_sans_ext(countryiso)
      
    }
    else{
      countryiso = countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c')
    }
    # Number of times we'll go through the loop
    show_modal_spinner(text=paste0("Creating Stations from ",input$ddata," with name ", countryiso )) # show the modal window
    database_creation(model=input$ddata, country = countryiso, clip_method = input$clipping, maxstations = input$maxstations, resol=inputshinnyresol(input$resol)  )
    remove_modal_spinner() # remove it when done
    states2 <- sf::st_read(paste(countryiso,"randomsample.shp",sep = "/"))
    leafletProxy("map") %>%  addCircles(data = states2, lng = ~x, lat = ~y, weight = 3,popup = ~cell, opacity = 0.5, group = "Step 1") %>% 
      # Layers control
      addLayersControl(
        overlayGroups = c("Boundaries","Step 1"),
        options = layersControlOptions(collapsed = FALSE)
      )
    
  })
  
  
  
  observeEvent(input$RunsA, {
    
    if(length(input$filemap) > 1 )
      
    {
      countryiso = countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c')
      req(input$filemap)
      shpdf <- input$filemap
      tempdirname <- dirname(shpdf$datapath[1])
      
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      
      countryiso = shpdf$name[grep(pattern = "*.shp$", shpdf$name)] 
      countryiso = tools::file_path_sans_ext(countryiso)
      
    }
    else{
      countryiso = countrycode(input$icountry_list, origin = 'country.name', destination = 'iso3c')
    }
    oldwd <- getwd()
    
    show_modal_spinner(text=paste0("Step 1 ",input$ddata," with name ", countryiso )) # show the modal window
    # show the modal window
    database_creation(model=input$ddata, country = countryiso, clip_method = input$clipping, maxstations = input$maxstations, resol=inputshinnyresol(input$resol)  )
    remove_modal_spinner()
    #get the number of stations
    stations <- read.csv(paste0(countryiso,"/BaseDatosEstaciones.csv", sep=""))
    #plot
    showNotification(paste("Stations:",length(unique(stations$id_station)), " stations detected!"))
    if(length(unique(stations$id_station)) < 14){
      stations_message = " At least 14 stations are needed to run this script, it may fail."
    } else{
      stations_message = " The number of stations is enough to run the script."
    }
    show_modal_spinner(text=paste0("Step 2 ",input$ddata," with name ", countryiso, ", Stations: ",length(unique(stations$id_station) ),".", stations_message )) # show the modal window
    output_2 <- step2_variable_calculation(country =  countryiso)
    remove_modal_spinner()
    
    show_modal_spinner(text=paste0("Step 3 ",input$ddata," with name ", countryiso )) # show the modal window
    output_3 <- regionalization(output_2$BaseDatosEstaciones,output_2$BaseRegistrosPr)
    remove_modal_spinner()
    
    show_modal_spinner(text=paste0("Exploratory Step ",input$ddata," with name ", countryiso )) # show the modal window
    exploratory(BaseDatosEstacionesClust = output_3$BaseDatosEstacionesClust, country = countryiso, VarInter="CumSumDec")
    remove_modal_spinner()
    
    show_modal_spinner(text=paste0("Step 4 ",input$ddata," with name ", countryiso )) # show the modal window
    output4 <- regional_frequency_analysis(output_3$ClustLevels,output_3$NombreClusters,output_3$VarInter, output_3$BaseDatosEstacionesClust, output_2$BaseRegistrosPr, output_2$z,output_3$Hc, countryiso)
    remove_modal_spinner()
    
    show_modal_spinner(text=paste0("Step 5 ",input$ddata," with name ", countryiso )) # show the modal window
    output5 <- period_estimation(output4$BaseSummaryStatistics, output4$BaseDatosEstacionesClust, output4$BaseBaseRegiones, output4$ResumeTable, output4$BaseMediaCompleta, output4$BaseProporCeros, output4$Basep0bias, output4$Baserfitdist, output4$Baserfitpara)
    remove_modal_spinner()
    
    show_modal_spinner(text=paste0("Creating rasters...",input$ddata," with name ", countryiso )) # show the modal window
    period_mapping(output4$ResumeTable, output5$BaseModelMapCor, countryiso, output_2$Boundaries)
    remove_modal_spinner()
    
    
    show_modal_spinner(text=paste0("Plotting... ", countryiso )) # show the modal window
    
    states2 <- sf::st_read("randomsample.shp")
    r <- raster('Mapas/RP1.tif')
    r05 <- raster('Mapas/RP05.tif')
    r06 <- raster('Mapas/RP06.tif')
    r07 <- raster('Mapas/RP07.tif')
    r08 <- raster('Mapas/RP08.tif')
    r09 <- raster('Mapas/RP09.tif')
    
    
    #rp legend colors
    minrp <- min(c(r@data@min,r05@data@min,r06@data@min,r07@data@min,r08@data@min,r09@data@min))
    maxrp <- max(c(r@data@max,r05@data@max,r06@data@max,r07@data@max,r08@data@max,r09@data@max))
    if(minrp < 0){
      minrp <- 0
    }
    rangodatos <- seq(minrp,maxrp,length.out=15)
    pal <- colorNumeric("Blues", rangodatos,
                        na.color = "transparent")
    
    
    
    EstQuant_1_in_2yr <- raster('Mapas/EstQuant_1_in_2yr.tif')
    EstQuant_1_in_5yr <- raster('Mapas/EstQuant_1_in_5yr.tif')
    EstQuant_1_in_10yr <- raster('Mapas/EstQuant_1_in_10yr.tif')
    EstQuant_1_in_20yr <- raster('Mapas/EstQuant_1_in_20yr.tif')
    EstQuant_1_in_30yr <- raster('Mapas/EstQuant_1_in_30yr.tif')
    EstQuant_1_in_40yr <- raster('Mapas/EstQuant_1_in_40yr.tif')
    EstQuant_1_in_50yr <- raster('Mapas/EstQuant_1_in_50yr.tif')
    EstQuant_1_in_60yr <- raster('Mapas/EstQuant_1_in_60yr.tif')
    EstQuant_1_in_70yr <- raster('Mapas/EstQuant_1_in_70yr.tif')
    EstQuant_1_in_80yr <- raster('Mapas/EstQuant_1_in_80yr.tif')
    EstQuant_1_in_90yr <- raster('Mapas/EstQuant_1_in_90yr.tif')
    EstQuant_1_in_100yr <- raster('Mapas/EstQuant_1_in_100yr.tif')
    
    
    #rp legend colors
    minrp2 <- min(c(EstQuant_1_in_2yr@data@min,EstQuant_1_in_5yr@data@min,EstQuant_1_in_10yr@data@min,EstQuant_1_in_20yr@data@min,EstQuant_1_in_30yr@data@min,EstQuant_1_in_40yr@data@min,EstQuant_1_in_50yr@data@min,EstQuant_1_in_60yr@data@min,EstQuant_1_in_70yr@data@min,EstQuant_1_in_80yr@data@min,EstQuant_1_in_90yr@data@min,EstQuant_1_in_100yr@data@min))
    maxrp2 <- max(c(EstQuant_1_in_2yr@data@max,EstQuant_1_in_5yr@data@max,EstQuant_1_in_10yr@data@max,EstQuant_1_in_20yr@data@max,EstQuant_1_in_30yr@data@max,EstQuant_1_in_40yr@data@max,EstQuant_1_in_50yr@data@max,EstQuant_1_in_60yr@data@max,EstQuant_1_in_70yr@data@max,EstQuant_1_in_80yr@data@max,EstQuant_1_in_90yr@data@max,EstQuant_1_in_100yr@data@max))
    if(minrp2 < 0){
      minrp2 <- 0
    }
    rangodatos2 <- seq(minrp2,maxrp2,length.out=5)
    pal2 <- colorNumeric("YlOrRd", rangodatos2,
                         na.color = "transparent")
    
    
    
    leafletProxy("map") %>%  addCircles(data = states2, lng = ~x, lat = ~y, weight = 3,popup = ~cell, opacity = 0.5, group = "Step 1") %>% 
      
      addRasterImage(r, colors = pal, opacity = 0.8, group = "RP1") %>%
      leaflet::addLegend(pal = pal, values = rangodatos,
                         title = "RP",   position = "bottomleft") %>%
      leaflet::addLegend(pal = pal2, values = rangodatos2,
                         title = "EstQuant",   position = "bottomright") %>%
      addRasterImage(r05, colors = pal, opacity = 0.8, group = "RP05") %>%
      addRasterImage(r06, colors = pal, opacity = 0.8, group = "RP06") %>%
      addRasterImage(r07, colors = pal, opacity = 0.8, group = "RP07") %>%
      addRasterImage(r08, colors = pal, opacity = 0.8, group = "RP08") %>%
      addRasterImage(r09, colors = pal, opacity = 0.8, group = "RP09") %>%
      
      addRasterImage(EstQuant_1_in_2yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_2yr") %>%
      addRasterImage(EstQuant_1_in_5yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_5yr") %>%
      addRasterImage(EstQuant_1_in_10yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_10yr") %>%
      addRasterImage(EstQuant_1_in_20yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_20yr") %>%
      addRasterImage(EstQuant_1_in_30yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_30yr") %>%
      addRasterImage(EstQuant_1_in_40yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_40yr") %>%
      addRasterImage(EstQuant_1_in_50yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_50yr") %>%
      addRasterImage(EstQuant_1_in_60yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_60yr") %>%
      addRasterImage(EstQuant_1_in_70yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_70yr") %>%
      addRasterImage(EstQuant_1_in_80yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_80yr") %>%
      addRasterImage(EstQuant_1_in_90yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_90yr") %>%
      addRasterImage(EstQuant_1_in_100yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_100yr") %>%
      
      
      addLayersControl(
        overlayGroups = c("Boundaries","Step 1", "RP05","RP06","RP07","RP08","RP09","RP1","EstQuant_1_in_2yr", "EstQuant_1_in_5yr", "EstQuant_1_in_10yr","EstQuant_1_in_20yr","EstQuant_1_in_30yr","EstQuant_1_in_40yr","EstQuant_1_in_50yr","EstQuant_1_in_60yr","EstQuant_1_in_70yr","EstQuant_1_in_80yr","EstQuant_1_in_90yr","EstQuant_1_in_100yr"  ),
        options = layersControlOptions(collapsed = TRUE), position = "topright"
      )  %>% 
      hideGroup(c("RP1","RP09","RP08","RP07","RP06","EstQuant_1_in_2yr", "EstQuant_1_in_5yr", "EstQuant_1_in_10yr","EstQuant_1_in_20yr","EstQuant_1_in_30yr","EstQuant_1_in_40yr","EstQuant_1_in_50yr","EstQuant_1_in_60yr","EstQuant_1_in_70yr","EstQuant_1_in_80yr","EstQuant_1_in_90yr","EstQuant_1_in_100yr"))
    
    setwd(oldwd)
    remove_modal_spinner()
    
    
    
    
    
    
    

    
    
  })

  
  
  
  observeEvent(input$RunsV, {
    
    if(length(input$filemap2) > 1 )
      
    {
      countryiso = countrycode(input$icountry_list2, origin = 'country.name', destination = 'iso3c')
      req(input$filemap2)
      shpdf <- input$filemap2
      tempdirname <- dirname(shpdf$datapath[1])
      
      # Rename files
      for (i in 1:nrow(shpdf)) {
        file.rename(
          shpdf$datapath[i],
          paste0(tempdirname, "/", shpdf$name[i])
        )
      }
      
      countryiso = shpdf$name[grep(pattern = "*.shp$", shpdf$name)] 
      countryiso = tools::file_path_sans_ext(countryiso)
      
    }
    else{
      countryiso = countrycode(input$icountry_list2, origin = 'country.name', destination = 'iso3c')
    }
    oldwd <- getwd()
    setwd(countryiso)
    show_modal_spinner(text=paste0("plotting... ", countryiso )) # show the modal window
    
    states2 <-  raster::shapefile("randomsample.shp")
    r <- raster('Mapas/RP1.tif')
    r05 <- raster('Mapas/RP05.tif')
    r06 <- raster('Mapas/RP06.tif')
    r07 <- raster('Mapas/RP07.tif')
    r08 <- raster('Mapas/RP08.tif')
    r09 <- raster('Mapas/RP09.tif')
    
    
    #rp legend colors
    minrp <- min(c(r@data@min,r05@data@min,r06@data@min,r07@data@min,r08@data@min,r09@data@min))
    maxrp <- max(c(r@data@max,r05@data@max,r06@data@max,r07@data@max,r08@data@max,r09@data@max))
    
    if(minrp < 0){
      minrp <- 0
    }
    
    rangodatos <- seq(minrp,maxrp,length.out=10)
    pal <- colorNumeric("Blues", rangodatos,
                        na.color = "transparent")
    
    
    
    EstQuant_1_in_2yr <- raster('Mapas/EstQuant_1_in_2yr.tif')
    EstQuant_1_in_5yr <- raster('Mapas/EstQuant_1_in_5yr.tif')
    EstQuant_1_in_10yr <- raster('Mapas/EstQuant_1_in_10yr.tif')
    EstQuant_1_in_20yr <- raster('Mapas/EstQuant_1_in_20yr.tif')
    EstQuant_1_in_30yr <- raster('Mapas/EstQuant_1_in_30yr.tif')
    EstQuant_1_in_40yr <- raster('Mapas/EstQuant_1_in_40yr.tif')
    EstQuant_1_in_50yr <- raster('Mapas/EstQuant_1_in_50yr.tif')
    EstQuant_1_in_60yr <- raster('Mapas/EstQuant_1_in_60yr.tif')
    EstQuant_1_in_70yr <- raster('Mapas/EstQuant_1_in_70yr.tif')
    EstQuant_1_in_80yr <- raster('Mapas/EstQuant_1_in_80yr.tif')
    EstQuant_1_in_90yr <- raster('Mapas/EstQuant_1_in_90yr.tif')
    EstQuant_1_in_100yr <- raster('Mapas/EstQuant_1_in_100yr.tif')
    
    #rp legend colors
    minrp2 <- min(c(EstQuant_1_in_2yr@data@min,EstQuant_1_in_5yr@data@min,EstQuant_1_in_10yr@data@min,EstQuant_1_in_20yr@data@min,EstQuant_1_in_30yr@data@min,EstQuant_1_in_40yr@data@min,EstQuant_1_in_50yr@data@min,EstQuant_1_in_60yr@data@min,EstQuant_1_in_70yr@data@min,EstQuant_1_in_80yr@data@min,EstQuant_1_in_90yr@data@min,EstQuant_1_in_100yr@data@min))
    maxrp2 <- max(c(EstQuant_1_in_2yr@data@max,EstQuant_1_in_5yr@data@max,EstQuant_1_in_10yr@data@max,EstQuant_1_in_20yr@data@max,EstQuant_1_in_30yr@data@max,EstQuant_1_in_40yr@data@max,EstQuant_1_in_50yr@data@max,EstQuant_1_in_60yr@data@max,EstQuant_1_in_70yr@data@max,EstQuant_1_in_80yr@data@max,EstQuant_1_in_90yr@data@max,EstQuant_1_in_100yr@data@max))
    if(minrp2 < 0){
      minrp2 <- 0
    }
    rangodatos2 <- seq(minrp2,maxrp2,length.out=5)
    
    pal2 <- colorNumeric("YlOrRd", rangodatos2,
                        na.color = "transparent")
    
    
    
    leafletProxy("map2") %>%  
      addCircles(data = states2, lng = ~x, lat = ~y, weight = 3,popup = ~cell, opacity = 0.5, group = "Step 1") %>% 
      addRasterImage(r, colors = pal, opacity = 0.8, group = "RP1") %>%
      leaflet::addLegend(pal = pal, values = rangodatos,
               title = "RP",   position = "bottomleft") %>%
      leaflet::addLegend(pal = pal2, values = rangodatos2,
                         title = "EstQuant",   position = "bottomright") %>%
      addRasterImage(r05, colors = pal, opacity = 0.8, group = "RP05") %>%
      addRasterImage(r06, colors = pal, opacity = 0.8, group = "RP06") %>%
      addRasterImage(r07, colors = pal, opacity = 0.8, group = "RP07") %>%
      addRasterImage(r08, colors = pal, opacity = 0.8, group = "RP08") %>%
      addRasterImage(r09, colors = pal, opacity = 0.8, group = "RP09") %>%
      
      addRasterImage(EstQuant_1_in_2yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_2yr") %>%
      addRasterImage(EstQuant_1_in_5yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_5yr") %>%
      addRasterImage(EstQuant_1_in_10yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_10yr") %>%
      addRasterImage(EstQuant_1_in_20yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_20yr") %>%
      addRasterImage(EstQuant_1_in_30yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_30yr") %>%
      addRasterImage(EstQuant_1_in_40yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_40yr") %>%
      addRasterImage(EstQuant_1_in_50yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_50yr") %>%
      addRasterImage(EstQuant_1_in_60yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_60yr") %>%
      addRasterImage(EstQuant_1_in_70yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_70yr") %>%
      addRasterImage(EstQuant_1_in_80yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_80yr") %>%
      addRasterImage(EstQuant_1_in_90yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_90yr") %>%
      addRasterImage(EstQuant_1_in_100yr, colors = pal2, opacity = 0.8, group = "EstQuant_1_in_100yr") %>%
      
      
      addLayersControl(
        overlayGroups = c("Boundaries","Step 1", "RP05","RP06","RP07","RP08","RP09","RP1","EstQuant_1_in_2yr", "EstQuant_1_in_5yr", "EstQuant_1_in_10yr","EstQuant_1_in_20yr","EstQuant_1_in_30yr","EstQuant_1_in_40yr","EstQuant_1_in_50yr","EstQuant_1_in_60yr","EstQuant_1_in_70yr","EstQuant_1_in_80yr","EstQuant_1_in_90yr","EstQuant_1_in_100yr"  ),
        options = layersControlOptions(collapsed = TRUE), position = "topright"
      )  %>% 
      hideGroup(c("RP1","RP09","RP08","RP07","RP06","EstQuant_1_in_2yr", "EstQuant_1_in_5yr", "EstQuant_1_in_10yr","EstQuant_1_in_20yr","EstQuant_1_in_30yr","EstQuant_1_in_40yr","EstQuant_1_in_50yr","EstQuant_1_in_60yr","EstQuant_1_in_70yr","EstQuant_1_in_80yr","EstQuant_1_in_90yr","EstQuant_1_in_100yr"))
    
    setwd(oldwd)
    remove_modal_spinner()
    
    
  })
  
  
}

# shinyApp()
shinyApp(ui = ui, server = server)