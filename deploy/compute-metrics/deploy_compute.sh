#!/bin/bash

#conda create --name prep_env python=3.9.7 --yes
#source activate prep_env
#conda install -c conda-forge pygeos=0.12.0 --yes
#conda install -c conda-forge geopandas=0.10.2 --yes
#conda install -c conda-forge urlpath=1.2.0 --yes
#conda install -c conda-forge pyarrow=6.0.1 --yes
#conda install -c conda-forge dask-geopandas=0.83 --yes

#country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)
country_chunk=(SLE)

log_file_arg=/Users/nm/Downloads/production/jobs/log_output.log
country_chunk_arg=${country_chunk[@]}
osm_dir_arg=/Users/nm/Downloads/production/inputs/osm
gadm_dir_arg=/Users/nm/Downloads/production/inputs/gadm
output_dir_arg=/Users/nm/Downloads/production/outputs

python batch_blocks.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg --output_dir $output_dir_arg
