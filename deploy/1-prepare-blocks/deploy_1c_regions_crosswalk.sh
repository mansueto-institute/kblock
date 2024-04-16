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

log_file_arg=$directory_path/jobs-reprex/deploy_1c_regions_crosswalk.log
country_chunk_arg=${iso_list[@]}
africapolis_file_arg=$directory_path/inputs-reprex/africapolis/AFRICAPOLIS2020.gpkg
ghsl_file_arg=$directory_path/inputs-reprex/ghsl/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg
gadm_dir_arg=$directory_path/inputs-reprex/gadm
blocks_dir_arg=$directory_path/outputs-reprex/blocks
output_dir_arg=$directory_path/outputs-reprex

python ./kblock/batch_1c_regions_crosswalk.py --log_file $log_file_arg --country_chunk $country_chunk_arg --africapolis_file $africapolis_file_arg --ghsl_file $ghsl_file_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg  --output_dir $output_dir_arg
