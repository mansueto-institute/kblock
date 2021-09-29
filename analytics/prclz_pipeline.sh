#!/bin/bash
#SBATCH --job-name=prclz_pipeline
#SBATCH --output=prclz_pipeline.ER
#SBATCH --error=prclz_pipeline.OU
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4000


## You'll need to change the directory names possibly, but the general idea is here


## To download the GADMS for every country, specify the repo, but use this general format:
## prclz download gadm {{your repo}}

## If you want geofabrik, run:

## prclz download geofabrik {{your repo}}

## To extract the geofabrik files, run:
## The geofabrik files are going to be the ones ending in .osm.pbf

for geofabrik_file in $(ls {{your repo}}; do
    geofabrik_path="{{your repo}}$geofabrik_file"
    prclz extract $geofabrik_path {{your chosen geofabrik extraction repo, could be the same as original geofabrik repo}}

## I'm not putting in anything about split buildings because I think Nico said that he's done that already, thanks Nico!

## To make the blocks, you'll need to match up the GADM files that end in _3 if that's the highest # or _4 if that's the 
## highest number, just choose the .shp GADM file in your GADM repo that has the highest _{{number}}. Match that up with
## the geofabrik lines geojson that got extracted in the last step. So, for SLE, it would be:

## prclz blocks gadm/SLE/gadm36_SLE_3.shp geofabrik/sierra-leone_lines.geojson blocks/

## It's not immediately obvious to me how to do that matching in Bash, so I'll leave it as pseudocode for now. If I were doing it in Python,
## If I were doing it in Python, I would create a CSV that matches every country code to its name, and loop through the rows of that
## csv.

## Loop through the buildings files and matched blocks files, calculating complexity. 
## Probably going to be better to create an array out of each, zip, and parallelize
## because TZA took > 24 hours. Or you could just run it multiple times without the 
## overwrite argument.

for buildings_file in $(ls /project2/bettencourt/mnp/prclz/data/buildings/Africa/TZA); do
    buildings_path="/project2/bettencourt/mnp/prclz/data/buildings/Africa/TZA/$buildings_file"
    blocks_file=${buildings_file##*buildings_}
    blocks_file=${blocks_file%.*}
    blocks_path="/project2/bettencourt/mnp/prclz/data/blocks/Africa/TZA/blocks_$blocks_file.csv"
    prclz complexity $blocks_path $buildings_path ~/prclz_test/complexity/
done

