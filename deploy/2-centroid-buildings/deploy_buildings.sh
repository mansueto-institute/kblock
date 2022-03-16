#!/bin/bash

#cd /Users/nm/Desktop/bldg_points
#source activate geo_zsh

country_list=(SLE LBR GIN)

log_file_arg=/Users/nm/Downloads/production/jobs/log_points.log
country_chunk_arg=${country_list[@]}
buildings_dir_arg=/Users/nm/Downloads/production/inputs/buildings/dev-dir
gadm_dir_arg=/Users/nm/Downloads/production/inputs/gadm
output_dir_arg=/Users/nm/Downloads/production/inputs/buildingpoints

python batch_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --buildings_dir $buildings_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg

