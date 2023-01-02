#!/bin/bash


#iso_list=(SLE GIN LBR DJI)

# iso_list=(DJI CPV COM STP SYC)

iso_list=(SYC DJI BDI)

working_directory=/Users/nm/Downloads/update

log_file_arg=$working_directory/jobs/deploy_1b_generate_blocks.log
country_chunk_arg=${iso_list[@]}
osm_dir_arg=$working_directory/inputs/osm/parquet
gadm_dir_arg=$working_directory/inputs/gadm/parquet
output_dir_arg=$working_directory/outputs

python batch_1b_generate_blocks.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg --output_dir $output_dir_arg

