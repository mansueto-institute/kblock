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

country=SLE

log_file_arg=/Users/MMSMITH/Desktop/logs/jobs/ghsl_log.log
ghsl_file=/Users/MMSMITH/Desktop/data/GHSL/buffered_GHSL.parquet
blocks_dir_arg=/Users/MMSMITH/Desktop/data/blocks
main_output_dir_arg=/Users/MMSMITH/Desktop/data/outputs
ancillary_output_dir_arg=/Users/MMSMITH/Desktop/data/outputs/ancillary

python batch_ghsl.py --log_file $log_file_arg --country $country --ghsl_file $ghsl_file --blocks_dir $blocks_dir_arg --output_dir $main_output_dir_arg --ancillary_dir $ancillary_output_dir_arg

