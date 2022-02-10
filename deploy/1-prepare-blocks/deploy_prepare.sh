#!/bin/bash

# conda create --name pygisparallel python=3.9.7 --yes
# source activate pygisparallel 
# conda install -c conda-forge pygeos=0.10.2 --yes
# conda install -c conda-forge geopandas=0.10.2 --yes 
# conda install -c conda-forge urlpath --yes
# conda install -c conda-forge dask --yes
# conda install -c conda-forge dask-geopandas --yes
# conda install -c conda-forge pyarrow --yes
# conda install -c conda-forge mpi4py --yes
# conda install -c conda-forge dask-mpi --yes
# conda install -c conda-forge dask-jobqueue --yes
# conda install -c conda-forge multiprocess --yes

country_chunk=(SLE DJI)

log_file_arg=/Users/nm/Downloads/production/jobs/log_output.log
country_chunk_arg=${country_chunk[@]}
osm_dir_arg=/Users/nm/Downloads/production/inputs/osm
gadm_dir_arg=/Users/nm/Downloads/production/inputs/gadm
output_dir_arg=/Users/nm/Downloads/production/outputs

python batch_prepare.py --log_file $log_file_arg --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg --output_dir $output_dir_arg
