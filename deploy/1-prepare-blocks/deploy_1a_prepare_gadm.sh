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

log_file_arg=$directory_path/jobs-reprex/deploy_1a_prepare_gadm.log
country_list_arg=${iso_list[@]}
gadm_dir_arg=$directory_path/inputs-reprex/gadm/geojson
daylight_dir_arg=$directory_path/inputs-reprex/daylight/coastlines-v1.19
osm_dir_arg=$directory_path/inputs-reprex/osm/parquet
output_dir_arg=$directory_path/inputs-reprex/gadm

python ./kblock/batch_1a_prepare_gadm.py --log_file $log_file_arg --country_chunk $country_list_arg --gadm_dir $gadm_dir_arg --daylight_dir $daylight_dir_arg --osm_dir $osm_dir_arg --output_dir $output_dir_arg

