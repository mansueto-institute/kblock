# etl

## __init__.py

Allows for the exportation of the functions under etl/

## country_codes.csv

A CSV that translates between country codes and geofabrik and gadm names so as to acquire the geofabrik and gadm URLS from the passed in country codes

## download.py

Contains functions that download GADM and geofabrik data. The main one, and the only one you should ever need to call, is called download. You can import it using '''from kblock.dev.etl.download import download'''. We can also add download to the __init__.py files in order to make it available from a higher level. 

### download

Call this function by passing the following: 

1. The string 'gadm' or 'geofabrik', depending on which you want to download

2. The directory where you want the downloaded data to be stored. Ideally a path from root for clarity. It seems to not work if you use a tilde to represent your home directory. 

3. A list of one or more country code strings representing the countries that you want to download the gadm or geofabrik files of. Can be any country code found in country_codes.csv. If you want to download multiple, specify them in a comma separated list without spaces. 

4. An optional boolean that controls whether your download will overwrite older versions of the same file. Defaults to False.

5. An optional boolean that controls whether your download will have a small progress bar. Default to False.

6. An optional boolean that controls whether you will be asked for confirmation as to whether you want to download to the directory you requested, with a warning that new directories may be created. Defaults to False. 


## extract.py

A single function script with the following function:

### extract

This function takes in a geofabrik pbf file and extracts its building polygons, building linestrings, and lines and saves them. Call this function by passing it: 

1. The file path of the pbf file. Ideally a path from root for clarity. It seems to not work if you use a tilde to represent your home directory. 

2. The folder you want the output files to be written to. 

3. An optional boolean that controls whether your extraction will overwrite older versions of the same file. Defaults to False. 


## extract.sh

This is the engine behind extract.py. It uses gdal's ogr2ogr to translate the pbf file into separate geojsons.

## osmconf.ini

This file is a necessary piece for ogr2ogr, it's basically a settings file for it. It must be in the same directory as extract.sh

