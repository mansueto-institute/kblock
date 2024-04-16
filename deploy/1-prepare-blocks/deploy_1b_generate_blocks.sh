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

log_file_arg=$directory_path/jobs-reprex/deploy_1b_generate_blocks.log
country_chunk_arg=${iso_list[@]}
osm_dir_arg=$directory_path/inputs-reprex/osm/parquet
gadm_dir_arg=$directory_path/inputs-reprex/gadm/parquet
output_dir_arg=$directory_path/outputs-reprex

python ./kblock/batch_1b_generate_blocks.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg --output_dir $output_dir_arg

