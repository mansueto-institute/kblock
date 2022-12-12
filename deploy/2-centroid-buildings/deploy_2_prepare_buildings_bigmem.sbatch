#!/bin/bash

#SBATCH --job-name=deploy_2_prepare_buildings_big
#SBATCH --partition=bigmem2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=320000
#SBATCH --output=/project2/bettencourt/mnp/update/jobs/deploy_2_prepare_buildings_bigmem.out
#SBATCH --error=/project2/bettencourt/mnp/update/jobs/deploy_2_prepare_buildings_bigmem.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nmarchio@uchicago.edu
#SBATCH --time=36:00:00
#SBATCH --account=pi-bettencourt

module load python/anaconda-2021.05
source activate geospatial

iso_list=(NGA ZAF SDN TZA UGA)

working_directory=/project2/bettencourt/mnp/update

codes_file_arg=$working_directory/repos/kblock/deploy/2-centroid-buildings/ecopia_country_codes.csv

log_file_arg=$working_directory/jobs/deploy_2_prepare_buildings_bigmem.log
progress_file_arg=$working_directory/jobs/ecopia_progress_bigmem.csv

country_chunk_arg=${iso_list[@]}
buildings_dir_arg=$working_directory/inputs/ecopia
blocks_dir_arg=$working_directory/outputs/blocks
output_dir_arg=$working_directory/outputs

python batch_2_prepare_buildings.py --log_file $log_file_arg  --country_chunk $country_chunk_arg --codes_file $codes_file_arg --progress_file $progress_file_arg --blocks_dir $blocks_dir_arg --buildings_dir $buildings_dir_arg --output_dir $output_dir_arg

