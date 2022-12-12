#!/bin/bash

# iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

# iso_list=(SLE GIN LBR DJI)

iso_list=(DJI CPV COM STP SYC)

working_directory=/Users/nm/Downloads/update

log_file_arg=$working_directory/jobs/deploy_1a_prepare_gadm.log
country_list_arg=${iso_list[@]}
gadm_dir_arg=$working_directory/inputs/gadm/geojson
daylight_dir_arg=$working_directory/inputs/daylight/coastlines-v1.19
osm_dir_arg=$working_directory/inputs/osm/parquet
output_dir_arg=$working_directory/inputs/gadm

python batch_1a_prepare_gadm.py --log_file $log_file_arg --country_chunk $country_list_arg --gadm_dir $gadm_dir_arg --daylight_dir $daylight_dir_arg --osm_dir $osm_dir_arg --output_dir $output_dir_arg

