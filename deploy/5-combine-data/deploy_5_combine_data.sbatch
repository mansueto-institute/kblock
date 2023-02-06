#!/bin/bash

#SBATCH --job-name=deploy_deploy_5_combine_data
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/deploy_5_combine_data.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/deploy_5_combine_data.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
source activate geospatial

iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER NGA RWA SDN SEN SLE SOM SSD STP SWZ SYC TCD TGO TZA UGA ZAF ZMB ZWE)

working_directory=/project2/bettencourt/mnp/update/outputs

log_file_arg=$working_directory/jobs/deploy_5_combine_data.log
country_chunk_arg=${iso_list[@]}
blocks_dir_arg=$working_directory/blocks
population_dir_arg=$working_directory/population
buildings_dir_arg=$working_directory/buildings
complexity_dir_arg=$working_directory/complexity
streets_dir_arg=$working_directory/streets
crosswalks_dir_arg=$working_directory/crosswalks
output_dir_arg=$working_directory/combined

python batch_5_combine_data.py --log_file $log_file_arg --country_chunk $country_chunk_arg --blocks_dir $blocks_dir_arg --population_dir $population_dir_arg --buildings_dir $buildings_dir_arg --complexity_dir $complexity_dir_arg --streets_dir $streets_dir_arg --crosswalks_dir $crosswalks_dir_arg --output_dir $output_dir_arg



 