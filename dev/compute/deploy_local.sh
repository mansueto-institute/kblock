#!/bin/bash

# conda create --name pygeos_env python=3.9.5
# source activate pygeos_env
# conda config --env --add channels conda-forge
# conda config --env --set channel_priority strict
# conda install pygeos --channel conda-forge
# conda install geopandas
# conda install gdal --channel conda-forge
# conda install -c conda-forge pyarrow
# conda install libpysal

gadm_list=(SLE.1.1.1_1 SLE.1.1.2_1 SLE.1.1.3_1 SLE.1.1.4_1 SLE.1.1.5_1 SLE.1.1.6_1)

log_file_arg=/Users/nm/Desktop/outputs/log_output.log
country_code_arg=SLE
country_code_file_arg=/Users/nm/Desktop/Projects/work/mnp/kblock/dev/compute/country_codes.csv
gadm_parent_dir_arg=/Users/nm/Desktop/Projects/work/mnp/mnp-dev/data/gadm
gadm_chunk_arg=${gadm_list[@]}
streets_parent_dir_arg=/Users/nm/Desktop/Projects/work/mnp/mnp-dev/data/osm
building_parent_dir_arg=/Users/nm/Downloads
output_dir_arg=/Users/nm/Desktop/outputs

python batch_kblock.py --log_file $log_file_arg --country_code $country_code_arg --country_code_file $country_code_file_arg --gadm_chunk $gadm_chunk_arg --gadm_parent_dir $gadm_parent_dir_arg --streets_parent_dir $streets_parent_dir_arg --building_parent_dir $building_parent_dir_arg --output_dir $output_dir_arg
