#!/bin/bash

#SBATCH --job-name=compute_metrics
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=54000
#SBATCH --nodes=3
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=14
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/compute_metrics_job.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/compute_metrics_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=24:00:00
#SBATCH --account=pi-bettencourt

module load openmpi/3.1.2 
module load python/anaconda-2021.05

source activate prep_env

mpirun python bath_compute.py log_file: Path, country_chunk: list, blocks_dir: Path, streets_dir: Path, buildings_dir: Path, raster_dir: Path, output_dir: Path)


#pip install mpi4py dask-mpi dask-ml


conda install -c conda-forge dask-mpi


country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA ZAF ZMB ZWE UGA RWA)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_build_blocks.log
country_chunk_arg=${country_list[@]}
osm_dir_arg=/project2/bettencourt/mnp/production/inputs/osm
gadm_dir_arg=/project2/bettencourt/mnp/production/inputs/gadm
output_dir_arg=/project2/bettencourt/mnp/production/outputs

python batch_prepare.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg
