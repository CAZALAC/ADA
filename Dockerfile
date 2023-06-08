FROM rocker/rstudio:4.3
#creating a working directory
RUN mkdir /srv/shiny-server/
#coping all needed files
COPY . /srv/shiny-server
RUN mkdir /home/test
COPY . /home/test
#COPY ./Extra/sources.list /home/analysis/sources.list
#updating sources list
#RUN sudo cp /srv/shiny-server/sources.list /etc/apt/sources.list
#Install Drought Atlas required packages 
RUN sudo apt-get update
RUN apt-get install g++ sqlite3 libsqlite3-dev libtiff5-dev curl pkg-config libjbig-dev proj-bin  gdal-bin -y  
RUN apt-get install libtiff5-dev libgdal-dev libproj-dev libexpat1-dev wx-common unixodbc-dev cmake -y
RUN apt-get install saga -y
#USER shiny
EXPOSE 3838
#RUN chmod -R 777 /srv/shiny-server/
#Install R packages
#RUN R -e "options(repos = \
# list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-01-06/')); \
#  TotalADAPackages=c('raster', 'rgdal', 'countrycode', 'rts', 'ncdf4','HelpersMG', 'latticeExtra', 'reshape2', 'gdata', 'gtools', 'plyr', 'sqldf', 'zoo', 'Kendall','zyp', 'car', 'rgeos', 'lmom' , 'lmomRFA','sp', 'rrcov', 'nsRFA', 'circular', 'reshape' , 'ModelMap', 'maptools', 'RSAGA' ,'stringr', 'rasterVis', 'hydroGOF' , 'randomForest','progress', 'deldir', 'ggplot2', 'edarf', 'chron','lattice', 'homtest'); \
#  install.packages(TotalADAPackages);"
#COPY InstallADApackages.R /home/analysis/InstallADApackages.R
#CMD R -e "source('/home/analysis/InstallADApackages.R')"
