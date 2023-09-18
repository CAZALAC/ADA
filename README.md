# African Drought Atlas

This repository contains the code and resources for the African Drought Atlas (ADA), a project aimed at developing a comprehensive understanding of hydro-climatic extreme events in the Africa region. The ADA is inspired by the successful Drought Atlas for Latin America and Caribbean released by UNESCO-IHP and CAZALAC.

## Contributing
We welcome contributions from researchers, meteorologists, hydrologists, and anyone interested in advancing our understanding of hydro-climatic extreme events in Africa. If you have valuable data, insights, or improvements to the methodology, we encourage you to contribute to this repository.

# Instructions to Set Up the ADA Project

These are the two available options to configure and run the ADA project. Choose the one that suits your needs best and enjoy working on the project!

## Option A (Recommended - Rstudio).

## Clone the Repository

First, clone the repository to your local machine using the following Git command:

```bash
git clone https://github.com/CAZALAC/ADA.git
```

### Download the Folder "ADAFolderPredictors"

Download the provided folder using the following link:

https://drive.google.com/drive/folders/1138qEcC7yOGRirtGFBT4WctSN1ef1mg9?usp=drive_link

Drag and drop the downloaded folder into the root script directory, or manually copy and paste it.

### Create a new project in RStudio.
Create a new project in RStudio and, in the RStudio console, execute the following command to restore the necessary packages:


```r
renv::restore()
```

Open `app.R` and Run the application using the `Run App` function.

Note: If you encounter any errors, you can try the following:

```r
renv::install()
```

## Option B (Docker - Rstudio-Server)

If you prefer to use Docker:

1. Ensure that Docker is installed on your system.

2. Execute the following command:

```bash
docker pull pabrojast/ada:0.2
```


Note: (OPTIONAL) If you want to build your own image, run the following command in the script root:

```bash
docker-compose up
```

If you have any further questions or need more assistance, please don't hesitate to ask! 

