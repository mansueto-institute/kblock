#!/bin/bash

#iso_list=(SLE LBR GIN DJI)

iso_list=(DJI CPV COM STP SYC)

working_directory=/Users/nm/Downloads/update

codes_file_arg=/Users/nm/Desktop/kblock/deploy/2-centroid-buildings/ecopia_country_codes.csv

log_file_arg=$working_directory/jobs/deploy_2_prepare_buildings.log
progress_file_arg=$working_directory/jobs/ecopia_progress.csv

country_chunk_arg=${iso_list[@]}
buildings_dir_arg=$working_directory/inputs/ecopia
blocks_dir_arg=$working_directory/outputs/blocks
output_dir_arg=$working_directory/outputs

python batch_2_prepare_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --codes_file $codes_file_arg --progress_file $progress_file_arg --blocks_dir $blocks_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg



