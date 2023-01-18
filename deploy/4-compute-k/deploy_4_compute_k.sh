#!/bin/bash

country_chunk=(SYC ) #DJI COM STP LBR CPV

working_directory=/Users/nm/Downloads/update

log_file_arg=$working_directory/jobs/deploy_4_compute_k.log
country_chunk_arg=${country_chunk[@]}
chunk_size_arg=10000
core_count_arg=1
blocks_dir_arg=$working_directory/outputs/blocks
streets_dir_arg=$working_directory/outputs/streets
buildings_dir_arg=$working_directory/outputs/buildings/points
output_dir_arg=$working_directory/outputs

python batch_4_compute_k.py --log_file $log_file_arg --country_chunk $country_chunk_arg --chunk_size $chunk_size_arg --core_count $core_count_arg --blocks_dir $blocks_dir_arg --streets_dir $streets_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

