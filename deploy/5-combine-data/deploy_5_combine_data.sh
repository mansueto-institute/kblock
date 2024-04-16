#!/bin/bash

iso_list=(DJI)

if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

directory_path=$1

if [ ! -d "$directory_path" ]; then
    echo "Error: '$directory_path' is not a directory."
    exit 1
fi

log_file_arg=$directory_path/jobs-reprex/deploy_5_combine_data.log
country_chunk_arg=${iso_list[@]}
blocks_dir_arg=$directory_path/outputs-reprex/blocks
population_dir_arg=$directory_path/outputs-reprex/population
buildings_dir_arg=$directory_path/outputs-reprex/buildings
complexity_dir_arg=$directory_path/outputs-reprex/complexity
streets_dir_arg=$directory_path/outputs-reprex/streets
crosswalks_dir_arg=$directory_path/outputs-reprex/crosswalks
output_dir_arg=$directory_path/outputs-reprex/combined

python ./kblock/batch_5_combine_data.py --log_file $log_file_arg --country_chunk $country_chunk_arg --blocks_dir $blocks_dir_arg --population_dir $population_dir_arg --buildings_dir $buildings_dir_arg --complexity_dir $complexity_dir_arg --streets_dir $streets_dir_arg --crosswalks_dir $crosswalks_dir_arg --output_dir $output_dir_arg



 