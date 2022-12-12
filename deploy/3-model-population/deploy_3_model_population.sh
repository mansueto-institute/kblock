#!/bin/bash

iso_list=(DJI CPV COM STP SYC)

working_directory=/Users/nm/Downloads/update

log_file_arg=$working_directory/jobs/deploy_3_model_population.log
country_chunk_arg=${iso_list[@]}
gadm_dir_arg=$working_directory/inputs/gadm/parquet
rasters_dir_arg=$working_directory/inputs/rasters
blocks_dir_arg=$working_directory/outputs/blocks
buildings_dir_arg=$working_directory/outputs/buildings/points
output_dir_arg=$working_directory/outputs

python batch_3_model_population.py --log_file $log_file_arg --country_chunk $country_chunk_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg --rasters_dir $rasters_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

