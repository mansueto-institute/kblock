#!/bin/bash

#SBATCH --job-name=deploy_2_prepare_buildings
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=56000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/deploy_2_prepare_buildings.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/deploy_2_prepare_buildings.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
source activate geospatial

iso_list=(AGO BDI BEN BFA BWA CAF CIV CMR COD COG COM CPV DJI ERI ESH ETH GAB GHA GIN GMB GNB GNQ KEN LBR LSO MDG MLI MOZ MRT MUS MWI NAM NER RWA SEN SLE SOM SSD STP SWZ SYC TCD TGO ZMB ZWE)

working_directory=/project2/bettencourt/mnp/update

codes_file_arg=$working_directory/repos/kblock/deploy/2-centroid-buildings/ecopia_country_codes.csv

log_file_arg=$working_directory/jobs/deploy_2_prepare_buildings.log
progress_file_arg=$working_directory/jobs/ecopia_progress.csv

country_chunk_arg=${iso_list[@]}
buildings_dir_arg=$working_directory/inputs/ecopia
blocks_dir_arg=$working_directory/outputs/blocks
output_dir_arg=$working_directory/outputs

python batch_2_prepare_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --codes_file $codes_file_arg --progress_file $progress_file_arg --blocks_dir $blocks_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg
