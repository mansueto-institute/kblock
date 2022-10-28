#!/bin/bash

# GHSL URL: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_STAT_UCDB2015MT_GLOBE_R2019A/V1-2/

ghsl_dir_arg=/Users/nm/Downloads/outputs/ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg
blocks_dir_arg=/Users/nm/Downloads/outputs/blocks
output_dir_arg=/Users/nm/Downloads/outputs/regions

python /Users/nm/Desktop/kblock/kblock/batch_regions.py --ghsl_dir $ghsl_dir_arg --blocks_dir $blocks_dir_arg --output_dir $output_dir_arg

