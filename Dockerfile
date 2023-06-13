FROM rocker/rstudio:4.3
#Install Drought Atlas required packages 
RUN sudo apt-get update && apt-get -y install --no-install-recommends \ 
g++ sqlite3 libsqlite3-dev libtiff5-dev curl pkg-config libjbig-dev proj-bin  gdal-bin \
libtiff5-dev libgdal-dev libproj-dev libexpat1-dev wx-common unixodbc-dev cmake \
saga libudunits2-dev 
#creating a working directory
RUN mkdir /home/rstudio/ADA
#coping all needed files
COPY . /home/rstudio/ADA
#Install R packages
RUN R -e "source('/home/rstudio/ADA/Extra/1.Install_packages.R')"
RUN chown -R rstudio /home/rstudio/ADA