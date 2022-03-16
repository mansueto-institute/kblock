#!/bin/bash

#SBATCH --job-name=popmodel
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/production/jobs/population_job.out
#SBATCH --error=/project2/bettencourt/mnp/production/jobs/population_job.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
#source activate pygisparallel

country_list=(SLE AGO NGA TZA ZAF ZMB BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER SDN SEN SOM SSD STP SWZ SYC TCD TGO ZWE UGA RWA)

log_file_arg=/project2/bettencourt/mnp/production/jobs/log_pop.log
country_chunk_arg=${country_list[@]}
gadm_dir_arg=/project2/bettencourt/mnp/production/inputs/gadm
blocks_dir_arg=/project2/bettencourt/mnp/production/outputs/blocks
population_dir_arg=/project2/bettencourt/mnp/production/inputs/population
buildings_dir_arg=/project2/bettencourt/mnp/production/inputs/buildingpoints
output_dir_arg=/project2/bettencourt/mnp/production/outputs/population

python batch_popmodel.py --log_file $log_file_arg --country_chunk $country_chunk_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg --population_dir $population_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

