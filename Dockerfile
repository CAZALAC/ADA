FROM rocker/rstudio:4.3
#creating a working directory
RUN mkdir /home/rstudio/ADA
#coping all needed files
COPY . /home/rstudio/ADA
#Install Drought Atlas required packages 
RUN sudo apt-get update
RUN apt-get install g++ sqlite3 libsqlite3-dev libtiff5-dev curl pkg-config libjbig-dev proj-bin  gdal-bin -y  
RUN apt-get install libtiff5-dev libgdal-dev libproj-dev libexpat1-dev wx-common unixodbc-dev cmake -y
RUN apt-get install saga -y
RUN R -e "source('/home/rstudio/ADA/Extra/1.Install_packages.R')"

#Install R packages