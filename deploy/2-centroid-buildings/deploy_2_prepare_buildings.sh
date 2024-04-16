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

codes_file_arg=$directory_path/inputs-reprex/ecopia/ecopia_country_codes.csv
log_file_arg=$directory_path/jobs-reprex/deploy_2_prepare_buildings.log
progress_file_arg=$directory_path/jobs-reprex/ecopia_progress.csv
country_chunk_arg=${iso_list[@]}
buildings_dir_arg=$directory_path/inputs-reprex/ecopia
blocks_dir_arg=$directory_path/outputs-reprex/blocks
output_dir_arg=$directory_path/outputs-reprex

python ./kblock/batch_2_prepare_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --codes_file $codes_file_arg --progress_file $progress_file_arg --blocks_dir $blocks_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg



