#!/bin/bash

country_list=(SLE COG)

log_file_arg=/Users/nm/Downloads/production/jobs/log_pop.log
country_chunk_arg=${country_list[@]}
gadm_dir_arg=/Users/nm/Downloads/production/inputs/gadm
blocks_dir_arg=/Users/nm/Downloads/production/outputs/blocks
population_dir_arg=/Users/nm/Downloads/production/inputs/population/tifs
buildings_dir_arg=/Users/nm/Downloads/production/inputs/buildingpoints
output_dir_arg=/Users/nm/Downloads/production/outputs/population

python batch_popmodel.py --log_file $log_file_arg --country_chunk $country_chunk_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg --population_dir $population_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

