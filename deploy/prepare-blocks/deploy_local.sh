#!/bin/bash

# conda create --name pygeospatial python=3.9.7 --yes
# source activate pygeospatial 
# conda install -c conda-forge pygeos=0.12.0 --yes
# conda install -c conda-forge geopandas=0.10.2 --yes 
# conda install -c conda-forge urlpath --yes
# conda install -c conda-forge dask --yes
# conda install -c conda-forge dask-geopandas --yes
# conda install -c conda-forge pyarrow --yes
# conda install -c conda-forge mpi4py --yes
# conda install -c conda-forge dask-mpi --yes
# conda install -c conda-forge dask-jobqueue --yes
# conda install -c conda-forge multiprocess --yes

#country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)
country_chunk=(SLE)

log_file_arg=/Users/nm/Downloads/production/jobs/log_output.log
country_chunk_arg=${country_chunk[@]}
osm_dir_arg=/Users/nm/Downloads/production/inputs/osm
gadm_dir_arg=/Users/nm/Downloads/production/inputs/gadm
output_dir_arg=/Users/nm/Downloads/production/outputs

mpirun -np 5 python batch_prepare.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg --output_dir $output_dir_arg

