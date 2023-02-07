#!/bin/bash


iso_list=(DJI SYC SLE)

working_directory=/Users/nm/Downloads/combine-dev

log_file_arg=/Users/nm/Downloads/update/jobs/deploy_5_combine_data.log
country_chunk_arg=${iso_list[@]}
blocks_dir_arg=$working_directory/blocks
population_dir_arg=$working_directory/population
buildings_dir_arg=$working_directory/buildings
complexity_dir_arg=$working_directory/complexity
streets_dir_arg=$working_directory/streets
crosswalks_dir_arg=$working_directory/crosswalks
output_dir_arg=$working_directory/combined

python batch_5_combine_data.py --log_file $log_file_arg --country_chunk $country_chunk_arg --blocks_dir $blocks_dir_arg --population_dir $population_dir_arg --buildings_dir $buildings_dir_arg --complexity_dir $complexity_dir_arg --streets_dir $streets_dir_arg --crosswalks_dir $crosswalks_dir_arg --output_dir $output_dir_arg



 