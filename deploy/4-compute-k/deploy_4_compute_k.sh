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

log_file_arg=$directory_path/jobs-reprex/deploy_4_compute_k.log
country_chunk_arg=${iso_list[@]}
chunk_size_arg=1000
core_count_arg=2
blocks_dir_arg=$directory_path/outputs-reprex/blocks
streets_dir_arg=$directory_path/outputs-reprex/streets
buildings_dir_arg=$directory_path/outputs-reprex/buildings/points
output_dir_arg=$directory_path/outputs-reprex

python ./kblock/batch_4_compute_k.py --log_file $log_file_arg --country_chunk $country_chunk_arg --chunk_size $chunk_size_arg --core_count $core_count_arg --blocks_dir $blocks_dir_arg --streets_dir $streets_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

