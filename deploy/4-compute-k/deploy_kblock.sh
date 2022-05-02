#!/bin/bash

# conda create --name geospatial python=3.9.7 --yes
# source activate geospatial
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
# conda install -c conda-forge rasterio --yes
# conda install -c conda-forge xarray --yes
# conda install -c conda-forge rioxarray --yes


country_chunk=(DJI) #COG COG ZAF DJI SYC STP LBR

log_file_arg=/Users/nm/Downloads/production/jobs/log_k.log
country_chunk_arg=${country_chunk[@]}
chunk_size_arg=200000
core_count_arg=3
blocks_dir_arg=/Users/nm/Downloads/production/outputs/blocks
streets_dir_arg=/Users/nm/Downloads/production/outputs/streets
buildings_dir_arg=/Users/nm/Downloads/production/inputs/buildingpoints
population_dir_arg=/Users/nm/Downloads/production/outputs/population
output_dir_arg=/Users/nm/Downloads/production/outputs

python batch_kblock.py --log_file $log_file_arg --country_chunk $country_chunk_arg --chunk_size $chunk_size_arg --core_count $core_count_arg --blocks_dir $blocks_dir_arg --streets_dir $streets_dir_arg --buildings_dir $buildings_dir_arg --population_dir $population_dir_arg --output_dir $output_dir_arg

