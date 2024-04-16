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

log_file_arg=$directory_path/jobs-reprex/deploy_3_model_population.log
country_chunk_arg=${iso_list[@]}
gadm_dir_arg=$directory_path/inputs-reprex/gadm/parquet
rasters_dir_arg=$directory_path/inputs-reprex/rasters
targets_file_arg=$directory_path/inputs-reprex/un/WPP2022_Demographic_Indicators_Medium.csv
blocks_dir_arg=$directory_path/outputs-reprex/blocks
buildings_dir_arg=$directory_path/outputs-reprex/buildings/points
output_dir_arg=$directory_path/outputs-reprex

python ./kblock/batch_3_model_population.py --log_file $log_file_arg --country_chunk $country_chunk_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg --rasters_dir $rasters_dir_arg --targets_file $targets_file_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

