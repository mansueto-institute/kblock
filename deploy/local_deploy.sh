#!/bin/bash

#conda create --name download_prep_env python=3.9.7 --yes
#source activate download_prep_env
#conda install -c conda-forge pygeos=0.12.0 --yes
#conda install -c conda-forge geopandas=0.10.2 --yes
#conda install -c conda-forge cykhash=1.0.2 --yes
#conda install -c conda-forge pyrosm=0.6.1 --yes
#conda install -c conda-forge urlpath --yes
#conda install -c conda-forge pyarrow --yes
#conda install -c conda-forge dask --yes
#conda install -c conda-forge dask-geopandas --yes

#country_list=(SYC SLE COM)

country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

log_file_arg=/Users/nm/Desktop/output/log_output.log
country_chunk_arg=${country_list[@]}
output_dir_arg=/Users/nm/Desktop/output

python batch_etl.py --log_file $log_file_arg --country_chunk $country_chunk_arg --output_dir $output_dir_arg
