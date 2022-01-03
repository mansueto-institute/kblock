#!/bin/bash

#SBATCH --job-name=build_block
#SBATCH --partition=build
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32000
#SBATCH --output=/project2/bettencourt/mnp/analytics/deployments/build_block_job.out
#SBATCH --error=/project2/bettencourt/mnp/analytics/deployments/build_block_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=8:00:00
#SBATCH --account=pi-bettencourt

country_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)
 
log_file_arg=/project2/bettencourt/mnp/analytics/deployments/log_build_blocks.log
country_chunk_arg=${country_list[@]}
output_dir_arg=/project2/bettencourt/mnp/analytics/outputs

python batch_etl.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --output_dir $output_dir_arg
