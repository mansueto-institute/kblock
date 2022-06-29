library(ggplot2)
library(sf)
library(tidyverse)
library(viridis)
library(patchwork)
library(scales)
library(tidygeocoder)
library(readxl)
library(units)
library(ggsn)
library(arrow)
library(sfarrow)
library(osmdata)

city_label = "Luanda"
country_label = 'Angola'
dir_path = './'
iso_code = 'AGO'
ghsl_delin = TRUE


# Grab link to zip file
ghsl_url <- 'http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.zip'

# Define future file path with new created folder (ghsl) that will be in temp directory
filedir <- paste0(tempdir(), '/ghsl/')
# Delete the file path
unlink(filedir, recursive = TRUE)
# Create the above directory/folder
dir.create(filedir)

# Define the URL as a character object
ghsl_dir <- paste0(ghsl_url)
# Download the zipped data finally
download.file(url = ghsl_url, destfile = paste0(filedir, basename(ghsl_dir)))
# Unzip the data
unzip(paste0(filedir,basename(ghsl_dir)), exdir= filedir)

city_data <- read_csv(list.files(path = paste0(filedir,
                                               list.files(path = filedir)[1]), 
                                 full.names = TRUE, all.files = TRUE)[6]) 


# Get UN pop data for Angola
location_url = 'https://population.un.org/wpp/Download/Files/4_Metadata/WPP2019_F01_LOCATIONS.XLSX'
tmp_filepath <- paste0(tempdir(), '/', basename(location_url))
download.file(url = paste0(location_url), destfile = tmp_filepath)
location_codes <- read_xlsx(path = tmp_filepath, sheet = 'Location', range = "A17:H306") 

location_codes <- location_codes %>%
  dplyr::select_all(~gsub("\\s+|\\.|\\/|,|\\*|-", "_", .)) %>%
  dplyr::rename_all(list(tolower)) %>%
  dplyr::filter(name == 'Country/Area') %>%
  dplyr::rename(country_name = region__subregion__country_or_area_) %>%
  dplyr::select(location_code, country_name, iso3_alpha_code)

un_population <- read_csv('https://population.un.org/wpp/Download/Files/1_Indicators%20(Standard)/CSV_FILES/WPP2019_TotalPopulationBySex.csv')






