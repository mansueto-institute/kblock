#!/bin/bash

#SBATCH --job-name=deploy_3_model_population
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/deploy_3_model_population.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/deploy_3_model_population.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
source activate geospatial

iso_list=(SLE AGO NGA TZA ZAF ZMB BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER SDN SEN SOM SSD STP SWZ SYC TCD TGO ZWE UGA RWA)

working_directory=/project2/bettencourt/mnp/update

log_file_arg=$working_directory/jobs/deploy_3_model_population.log
country_chunk_arg=${iso_list[@]}
gadm_dir_arg=$working_directory/inputs/gadm/parquet
rasters_dir_arg=$working_directory/inputs/rasters
blocks_dir_arg=$working_directory/outputs/blocks
buildings_dir_arg=$working_directory/outputs/buildings/points
output_dir_arg=$working_directory/outputs

python batch_3_model_population.py --log_file $log_file_arg --country_chunk $country_chunk_arg --gadm_dir $gadm_dir_arg --blocks_dir $blocks_dir_arg --rasters_dir $rasters_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg
