#!/bin/bash

# GHSL URL: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/

ghsl_dir_arg=/Users/nm/Downloads/production/inputs/ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg
blocks_dir_arg=/Users/nm/Downloads/outputs/blocks
output_dir_arg=/Users/nm/Downloads/production/outputs/crosswalk

python /Users/nm/Desktop/repos/kblock/kblock/batch_crosswalk.py --ghsl_dir $ghsl_dir_arg --blocks_dir $blocks_dir_arg  --output_dir $output_dir_arg

