#!/bin/bash

#SBATCH --job-name=build_block
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/build_block_job.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/build_block_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05

source activate pygisparallel

country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA ZAF ZMB ZWE UGA RWA)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_build_blocks.log
country_chunk_arg=${country_list[@]}
osm_dir_arg=/project2/bettencourt/mnp/production/inputs/osm
gadm_dir_arg=/project2/bettencourt/mnp/production/inputs/gadm
output_dir_arg=/project2/bettencourt/mnp/production/outputs

python batch_prepare.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --osm_dir $osm_dir_arg --gadm_dir $gadm_dir_arg  --output_dir $output_dir_arg

